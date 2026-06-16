#include "CalibCLI.hpp"

#include <algorithm>
#include <cmath>
#include <clocale>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/log/trivial.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "libslic3r/Flow.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/Utils.hpp"

namespace Slic3r {

CLICalibType cli_calib_type_from_string(const std::string &name)
{
    if (name == "flow-yolo-recommended")    return CLICalibType::FlowRate_YOLO_Recommended;
    if (name == "flow-yolo-perfectionist")  return CLICalibType::FlowRate_YOLO_Perfectionist;
    if (name == "flow-pass1")               return CLICalibType::FlowRate_Pass1;
    if (name == "flow-pass2")               return CLICalibType::FlowRate_Pass2;
    if (name == "z-offset-pattern")         return CLICalibType::ZOffsetPattern;
    return CLICalibType::NoCalib;
}

std::string cli_calib_resource_path(CLICalibType type)
{
    const std::string root = resources_dir();
    switch (type) {
    case CLICalibType::FlowRate_YOLO_Recommended:
        return root + "/calib/filament_flow/Orca-LinearFlow.3mf";
    case CLICalibType::FlowRate_YOLO_Perfectionist:
        return root + "/calib/filament_flow/Orca-LinearFlow_fine.3mf";
    case CLICalibType::FlowRate_Pass1:
        return root + "/calib/filament_flow/flowrate-test-pass1.3mf";
    case CLICalibType::FlowRate_Pass2:
        return root + "/calib/filament_flow/flowrate-test-pass2.3mf";
    case CLICalibType::ZOffsetPattern:
        // Use any existing 3MF as a load-pipeline placeholder; cli_build_zcal_pattern wipes
        // the Model and rebuilds from procedural primitives.
        return root + "/calib/filament_flow/Orca-LinearFlow.3mf";
    default:
        return std::string();
    }
}

void cli_flowrate_params_for_type(CLICalibType type, int &pass, bool &is_linear)
{
    switch (type) {
    case CLICalibType::FlowRate_YOLO_Recommended:   pass = 1; is_linear = true;  break;
    case CLICalibType::FlowRate_YOLO_Perfectionist: pass = 2; is_linear = true;  break;
    case CLICalibType::FlowRate_Pass1:              pass = 1; is_linear = false; break;
    case CLICalibType::FlowRate_Pass2:              pass = 2; is_linear = false; break;
    default:                                        pass = 1; is_linear = true;  break;
    }
}

InfillPattern cli_parse_flow_pattern(const std::string &name)
{
    if (name == "monotonic-line" || name == "monotonic" || name == "monotonicline")
        return ipMonotonicLine;
    // default — matches the GUI dialog's default selection
    return ipArchimedeanChords;
}

//ORCA: parse the modifier from a flow-rate calibration object name. Filename format:
//      "flowrate_0" / "flowrate_0.035" / "flowrate_m0.04" (m = minus). Returns the modifier
//      as a float; defaults to 1.0 on parse failure (matches GUI fallback at Plater.cpp:13054).
//      Sets ok=false on parse failure so callers (like the block filter) can skip the object.
static float parse_flow_modifier(const std::string &object_name, bool &ok)
{
    ok = false;
    if (object_name.length() <= 9) return 1.0f;
    std::string mod_str = object_name.substr(9);
    if (!mod_str.empty() && mod_str[0] == 'm')
        mod_str[0] = '-';
    const std::string saved_locale = std::setlocale(LC_NUMERIC, nullptr);
    std::setlocale(LC_NUMERIC, "C");
    float modifier = 1.0f;
    try {
        modifier = std::stof(mod_str);
        ok = true;
    } catch (...) {}
    std::setlocale(LC_NUMERIC, saved_locale.c_str());
    return modifier;
}

//ORCA: filter flow-rate calibration blocks by count and/or |modifier|. Removes objects whose
//      modifier is too far from zero. Both filters compose (intersection). Symmetric removal
//      keeps the cluster centered on the existing centroid, so callers don't need to re-center.
static void filter_flowrate_blocks(Model &model, int max_blocks, double max_modifier)
{
    const bool count_filter = (max_blocks > 0 && static_cast<int>(model.objects.size()) > max_blocks);
    const bool range_filter = (max_modifier > 0.0);
    if (!count_filter && !range_filter) return;

    struct Entry { size_t idx; float mod; float abs_mod; };
    std::vector<Entry> entries;
    entries.reserve(model.objects.size());
    for (size_t i = 0; i < model.objects.size(); ++i) {
        bool ok = false;
        const float mod = parse_flow_modifier(model.objects[i]->name, ok);
        if (!ok) continue; // unparseable name — leave in place (treat as keeper)
        entries.push_back({i, mod, std::fabs(mod)});
    }
    if (entries.empty()) return;

    std::vector<bool> keep(model.objects.size(), true);
    // Range filter: drop blocks where |modifier| > max_modifier.
    if (range_filter) {
        for (const auto& e : entries)
            if (e.abs_mod > max_modifier + 1e-9) keep[e.idx] = false;
    }
    // Count filter: among still-kept entries, keep only the N closest to 0.
    if (count_filter) {
        std::vector<Entry> kept_entries;
        for (const auto& e : entries)
            if (keep[e.idx]) kept_entries.push_back(e);
        std::sort(kept_entries.begin(), kept_entries.end(),
                  [](const Entry& a, const Entry& b) { return a.abs_mod < b.abs_mod; });
        for (size_t k = static_cast<size_t>(max_blocks); k < kept_entries.size(); ++k)
            keep[kept_entries[k].idx] = false;
    }

    // Erase from highest index down via Model::delete_object so destruction goes through the
    // friend channel (ModelObject's dtor is private to the Model factory). Indices stay valid
    // because we iterate top-down.
    size_t removed = 0;
    for (size_t i = model.objects.size(); i-- > 0; ) {
        if (!keep[i]) {
            model.delete_object(i);
            ++removed;
        }
    }
    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_flowrate_calib: block filter removed %1% objects (max_blocks=%2% max_modifier=%3$.4f), %4% remain")
        % removed % max_blocks % max_modifier % model.objects.size();
}

//ORCA: Headless port of Plater::adjust_settings_for_flowrate_calib (Plater.cpp:12967). The original
//      function scaled the model via wxGetApp().plater()->canvas3D()->get_selection().scale(...),
//      which is GUI-bound; here we apply the same scale factors directly via ModelObject::scale().
//      Every per-object config setter is line-for-line identical to the GUI implementation so the
//      sliced G-code is byte-equivalent to a GUI wizard run with the same parameters.
void cli_apply_flowrate_calib(Model &model, DynamicPrintConfig &full_config, const CLIFlowRateParams &params)
{
    if (params.pass != 1 && params.pass != 2) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_flowrate_calib: invalid pass " << params.pass;
        return;
    }

    //ORCA: subset the loaded block-pair objects per --flow-blocks / --flow-range before any
    //      per-object/per-config work. Symmetric (keeps blocks closest to modifier=0).
    filter_flowrate_blocks(model, params.max_blocks, params.max_modifier);

    const ConfigOptionFloats *nozzle_diameter_config = full_config.option<ConfigOptionFloats>("nozzle_diameter");
    if (!nozzle_diameter_config || nozzle_diameter_config->values.empty()) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_flowrate_calib: no nozzle_diameter in config";
        return;
    }
    const double nozzle_diameter = nozzle_diameter_config->values[0];
    const double xyScale         = nozzle_diameter / 0.6;
    const double layer_height    = nozzle_diameter / 2.0;
    double       first_layer_h   = full_config.option<ConfigOptionFloat>("initial_layer_print_height")->value;
    first_layer_h                = std::max(first_layer_h, layer_height);
    const double zscale          = (first_layer_h + 9 * layer_height) / 2;

    //ORCA: apply scale via direct ModelObject transformation (matches GUI's selection.scale() semantics).
    const Vec3d scale_v = (xyScale > 1.2) ? Vec3d(xyScale, xyScale, zscale) : Vec3d(1.0, 1.0, zscale);
    for (ModelObject *mo : model.objects) {
        if (!mo) continue;
        mo->scale(scale_v);
    }

    //ORCA: derive max_infill_speed exactly as the GUI does. cur_flowrate (filament_flow_ratio) and
    //      filament_max_volumetric_speed come from the filament preset slot 0.
    const double cur_flowrate              = full_config.option<ConfigOptionFloats>("filament_flow_ratio")->get_at(0);
    const Flow   infill_flow               = Flow(nozzle_diameter * 1.2f, layer_height, nozzle_diameter);
    const double filament_max_vol_speed    = full_config.option<ConfigOptionFloats>("filament_max_volumetric_speed")->get_at(0);
    double       max_infill_speed;
    if (params.is_linear) {
        max_infill_speed = filament_max_vol_speed /
                           (infill_flow.mm3_per_mm() * (cur_flowrate + (params.pass == 2 ? 0.035 : 0.05)) / cur_flowrate);
    } else {
        max_infill_speed = filament_max_vol_speed / (infill_flow.mm3_per_mm() * (params.pass == 1 ? 1.2 : 1));
    }
    const double internal_solid_speed = std::floor(std::min(full_config.opt_float("internal_solid_infill_speed"), max_infill_speed));
    const double top_surface_speed    = std::floor(std::min(full_config.opt_float("top_surface_speed"), max_infill_speed));

    //ORCA: per-object setters — line-for-line port of the GUI implementation. Order preserved.
    for (ModelObject *mo : model.objects) {
        if (!mo) continue;
        mo->ensure_on_bed();
        mo->config.set_key_value("wall_loops", new ConfigOptionInt(1));
        mo->config.set_key_value("only_one_wall_top", new ConfigOptionBool(true));
        mo->config.set_key_value("thick_internal_bridges", new ConfigOptionBool(false));
        mo->config.set_key_value("enable_extra_bridge_layer", new ConfigOptionEnum<EnableExtraBridgeLayer>(eblDisabled));
        mo->config.set_key_value("internal_bridge_density", new ConfigOptionPercent(100));
        mo->config.set_key_value("sparse_infill_density", new ConfigOptionPercent(35));
        mo->config.set_key_value("min_width_top_surface", new ConfigOptionFloatOrPercent(100, true));
        mo->config.set_key_value("bottom_shell_layers", new ConfigOptionInt(2));
        mo->config.set_key_value("top_shell_layers", new ConfigOptionInt(5));
        mo->config.set_key_value("top_shell_thickness", new ConfigOptionFloat(0));
        mo->config.set_key_value("bottom_shell_thickness", new ConfigOptionFloat(0));
        mo->config.set_key_value("detect_thin_wall", new ConfigOptionBool(true));
        mo->config.set_key_value("filter_out_gap_fill", new ConfigOptionFloat(0));
        mo->config.set_key_value("sparse_infill_pattern", new ConfigOptionEnum<InfillPattern>(ipRectilinear));
        mo->config.set_key_value("top_surface_line_width", new ConfigOptionFloatOrPercent(nozzle_diameter * 1.2f, false));
        mo->config.set_key_value("internal_solid_infill_line_width", new ConfigOptionFloatOrPercent(nozzle_diameter * 1.2f, false));
        mo->config.set_key_value("top_surface_pattern", new ConfigOptionEnum<InfillPattern>(params.pattern));
        mo->config.set_key_value("top_solid_infill_flow_ratio", new ConfigOptionFloat(1.0f));
        mo->config.set_key_value("infill_direction", new ConfigOptionFloat(45));
        mo->config.set_key_value("solid_infill_direction", new ConfigOptionFloat(135));
        mo->config.set_key_value("align_infill_direction_to_model", new ConfigOptionBool(true));
        mo->config.set_key_value("ironing_type", new ConfigOptionEnum<IroningType>(IroningType::NoIroning));
        mo->config.set_key_value("internal_solid_infill_speed", new ConfigOptionFloat(internal_solid_speed));
        mo->config.set_key_value("top_surface_speed", new ConfigOptionFloat(top_surface_speed));
        mo->config.set_key_value("seam_slope_type", new ConfigOptionEnum<SeamScarfType>(SeamScarfType::None));
        mo->config.set_key_value("gap_fill_target", new ConfigOptionEnum<GapFillTarget>(GapFillTarget::gftNowhere));
        mo->config.set_key_value("calib_flowrate_topinfill_special_order", new ConfigOptionBool(true));

        //ORCA: parse the per-block flow rate modifier from the object name (filename format
        //      "flowrate_xxx"). 'm' prefix means negative (e.g. "flowrate_m0.01" → -0.01).
        //      print_flow_ratio differs between linear (YOLO) and non-linear (Pass1/Pass2):
        //        linear:     (cur_flowrate + modifier) / cur_flowrate
        //        non-linear: 1.0 + modifier / 100
        //      Identical to the GUI logic at Plater.cpp:13044-13065.
        bool mod_ok = false;
        const float modifier = parse_flow_modifier(mo->name, mod_ok);
        if (mod_ok) {
            const float pfr = params.is_linear
                ? (static_cast<float>(cur_flowrate) + modifier) / static_cast<float>(cur_flowrate)
                : 1.0f + modifier / 100.0f;
            mo->config.set_key_value("print_flow_ratio", new ConfigOptionFloat(pfr));
        }

        //ORCA: brim toggle from PR #13548 — when enabled, switch each block-pair object to an outer
        //      brim with the requested width and zero object-gap (matching the dialog's behavior).
        if (params.brim_enabled) {
            mo->config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btOuterOnly));
            mo->config.set_key_value("brim_width", new ConfigOptionFloat(params.brim_width));
            mo->config.set_key_value("brim_object_gap", new ConfigOptionFloat(0.0));
        }
    }

    //ORCA: brim auto-spacing — line-for-line port of the GUI block at Plater.cpp:13076-13115.
    //      Native 3MF spacing is 31mm with 30mm blocks (1mm edge gap); a 2mm brim would merge
    //      adjacent block brims. Scale object positions outward from the cluster centroid by the
    //      smallest uniform factor that widens every gap to required_gap = 2*brim_width + extra_gap.
    if (params.brim_enabled && model.objects.size() >= 2) {
        std::vector<Vec2d>         centers;
        std::vector<BoundingBoxf3> bbs;
        centers.reserve(model.objects.size());
        bbs.reserve(model.objects.size());
        for (auto mo : model.objects) {
            if (!mo) continue;
            const BoundingBoxf3 bb = mo->instance_bounding_box(0);
            bbs.push_back(bb);
            centers.emplace_back(0.5 * (bb.min.x() + bb.max.x()), 0.5 * (bb.min.y() + bb.max.y()));
        }
        Vec2d centroid(0.0, 0.0);
        for (const auto& c : centers) centroid += c;
        centroid /= double(centers.size());

        double min_center_dist = std::numeric_limits<double>::infinity();
        double min_edge_gap    = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < model.objects.size(); ++i) {
            for (size_t j = i + 1; j < model.objects.size(); ++j) {
                const double cd = (centers[i] - centers[j]).norm();
                if (cd < min_center_dist) min_center_dist = cd;
                const double dx = std::max(0.0, std::max(bbs[i].min.x() - bbs[j].max.x(), bbs[j].min.x() - bbs[i].max.x()));
                const double dy = std::max(0.0, std::max(bbs[i].min.y() - bbs[j].max.y(), bbs[j].min.y() - bbs[i].max.y()));
                const double eg = std::max(dx, dy); // axis-aligned grid: the larger axial gap is the actual gap
                if (eg < min_edge_gap) min_edge_gap = eg;
            }
        }

        const double required_gap = 2.0 * params.brim_width + params.brim_extra_gap;
        if (min_center_dist > 0.0 && min_edge_gap < required_gap) {
            const double extra_needed = required_gap - min_edge_gap;
            const double scale        = 1.0 + extra_needed / min_center_dist;
            for (size_t i = 0; i < model.objects.size(); ++i) {
                const Vec2d new_center = centroid + scale * (centers[i] - centroid);
                const Vec2d delta      = new_center - centers[i];
                model.objects[i]->translate_instances(Vec3d(delta.x(), delta.y(), 0.0));
            }
            BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_flowrate_calib: brim auto-spacing required_gap=%1$.2f min_edge_gap=%2$.2f scale=%3$.4f")
                % required_gap % min_edge_gap % scale;
        }
    }

    //ORCA: print-config overrides — the GUI applies these at the END of
    //      adjust_settings_for_flowrate_calib (Plater.cpp:13117-13121, plus
    //      :13042 max_volumetric_extrusion_rate_slope inside the per-object loop).
    //      These must run after the per-object loop so they're not silently shadowed.
    full_config.set_key_value("layer_height", new ConfigOptionFloat(layer_height));
    full_config.set_key_value("initial_layer_print_height", new ConfigOptionFloat(first_layer_h));
    full_config.set_key_value("alternate_extra_wall", new ConfigOptionBool(false));
    full_config.set_key_value("reduce_crossing_wall", new ConfigOptionBool(true));
    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));
    full_config.set_key_value("max_volumetric_extrusion_rate_slope", new ConfigOptionFloat(0));
    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_flowrate_calib: pass=%1% is_linear=%2% pattern=%3% brim=%4% objects=%5%")
        % params.pass % params.is_linear % static_cast<int>(params.pattern) % params.brim_enabled % model.objects.size();
}

// ============================================================
// Z-offset / Live-Z calibration pattern (procedural mesh).
// ============================================================

namespace {

// Layer thickness for every primitive in the Z-cal plate. Matches initial_layer_print_height
// so each mesh slices to exactly one layer regardless of layer_height.
constexpr double ZCAL_LAYER_THICKNESS = 0.20;

// Line widths chosen so single-perimeter walls render as one extrusion at 0.4mm nozzle.
constexpr double ZCAL_THIN_WALL_WIDTH = 0.45;

// Centered axis-aligned box ModelObject. The mesh is built at the origin (size centered in
// XY, sitting on Z=0 with its top at thickness) and the instance offset positions it on the
// plate. All Z-cal primitives use this builder.
ModelObject* add_box_object(Model& model, const std::string& name,
                            double xc, double yc,
                            double w, double h, double thickness)
{
    // its_make_cube draws from (0,0,0) to (w,h,t); shift by -w/2,-h/2,0 to center XY around origin.
    TriangleMesh mesh(its_make_cube(w, h, thickness));
    mesh.translate(-w/2.0f, -h/2.0f, 0.0f);
    ModelObject* mo = model.add_object();
    mo->name = name;
    // add_volume with modify_to_center_geometry=false keeps our mesh centered as built — the
    // default would recompute and apply origin_translation which we don't want for a pre-built
    // primitive at a known position.
    mo->add_volume(std::move(mesh), ModelVolumeType::MODEL_PART, /*modify_to_center_geometry=*/false);
    ModelInstance* inst = mo->add_instance();
    inst->set_offset(Vec3d(xc, yc, 0.0));
    return mo;
}

// Per-object overrides that make a small thin box render as a single-perimeter freestanding wall
// (no infill, no top/bottom). Used for Zone W walls and Zone G gap-pair walls.
void set_single_wall_config(ModelObject* mo)
{
    mo->config.set_key_value("wall_loops", new ConfigOptionInt(1));
    mo->config.set_key_value("top_shell_layers", new ConfigOptionInt(0));
    mo->config.set_key_value("bottom_shell_layers", new ConfigOptionInt(1));
    mo->config.set_key_value("sparse_infill_density", new ConfigOptionPercent(0));
    mo->config.set_key_value("only_one_wall_top", new ConfigOptionBool(false));
}

// Per-object overrides for Zone S — solid pad with concentric infill, the squish-uniformity
// indicator. One wall, concentric top pattern, no sparse infill (single layer is all-top).
void set_solid_pad_config(ModelObject* mo)
{
    mo->config.set_key_value("wall_loops", new ConfigOptionInt(1));
    mo->config.set_key_value("top_shell_layers", new ConfigOptionInt(1));
    mo->config.set_key_value("bottom_shell_layers", new ConfigOptionInt(1));
    mo->config.set_key_value("sparse_infill_density", new ConfigOptionPercent(0));
    mo->config.set_key_value("top_surface_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentric));
    mo->config.set_key_value("bottom_surface_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentric));
}

// Per-object overrides for Zone C corner loops — small ring (single-perimeter) at each corner.
// Same as single_wall_config but kept as its own helper so future tuning (e.g. brim toggle for
// corner adhesion) doesn't drag Zone W along.
void set_corner_loop_config(ModelObject* mo)
{
    set_single_wall_config(mo);
}

} // anonymous

void cli_build_zcal_pattern(Model &model, DynamicPrintConfig &full_config, const CLIZCalParams &params)
{
    //ORCA: the load pipeline populated model from a placeholder 3MF; wipe everything and rebuild.
    while (!model.objects.empty())
        model.delete_object(model.objects.size() - 1);

    const double plate     = params.plate_size;       // typically 100mm
    const double zone      = params.zone_size;        // typically 30mm (or smaller if zones don't fit)
    const double half      = plate / 2.0;

    // Layout in model space (centered at 0,0). The CLI's center-on-bed step later translates
    // the whole plate to the bed center.
    //
    //   fid_TL                                           fid_TR
    //     +                                                 +
    //     C_TL                                           C_TR
    //
    //              [ S ]      [ G ]      [ W ]
    //
    //     C_BL                                           C_BR
    //     +     |||||||||||||||||| scale bar               +
    //   fid_BL                                           fid_BR
    //
    // Fiducials at ±(half-5), Zone-C loops at ±(half-12), main row of S/G/W at y=+8.

    const double fid_off       = half - 5.0;     // fiducial center inset from plate edge
    const double corner_off    = half - 12.0;    // Zone C loop inset
    const double row_y         = 8.0;            // main S/G/W row y center
    const double zone_spacing  = 3.0;            // gap between adjacent main zones
    const double zone_pitch    = zone + zone_spacing;
    // Three centers: -pitch, 0, +pitch
    const double zone_S_xc     = -zone_pitch;
    const double zone_G_xc     =  0.0;
    const double zone_W_xc     = +zone_pitch;
    const double scale_bar_y   = -(half - 8.0);  // scale bar along bottom inset 8mm from edge
    const double scale_bar_len = 10.0;            // 10mm scale
    const double scale_bar_x0  = -scale_bar_len / 2.0;

    // ---- Fiducials (4× filled square at known XY) ----
    //      Single mesh per fiducial to avoid gcode-path conflicts from overlapping ModelObjects
    //      (cross shape needed two perpendicular bars overlapping at the center). A 5×5 solid
    //      square is just as detectable by AI-vision corner-finders as a "+" — the alignment
    //      pipeline locates the centroid of each filled blob.
    if (params.fiducials) {
        struct { double x, y; const char* name; } fids[4] = {
            { -fid_off, -fid_off, "fid_BL" },
            { +fid_off, -fid_off, "fid_BR" },
            { +fid_off, +fid_off, "fid_TR" },
            { -fid_off, +fid_off, "fid_TL" },
        };
        const double fid_size = 5.0;
        for (auto &f : fids) {
            ModelObject* mo = add_box_object(model, f.name,
                                             f.x, f.y, fid_size, fid_size, ZCAL_LAYER_THICKNESS);
            set_solid_pad_config(mo);
        }
    }

    // ---- Scale bar (1 baseline + 11 ticks at 1mm intervals) ----
    if (params.scale_bar) {
        // Baseline: 10mm long, 0.5mm wide, at y=scale_bar_y.
        ModelObject* baseline = add_box_object(model, "scale_baseline",
                                                scale_bar_x0 + scale_bar_len / 2.0, scale_bar_y,
                                                scale_bar_len, 0.5, ZCAL_LAYER_THICKNESS);
        set_solid_pad_config(baseline);

        // 11 ticks: short verticals at x = scale_bar_x0 + 0, 1, 2, ..., 10mm.
        for (int i = 0; i <= 10; ++i) {
            const double tx = scale_bar_x0 + static_cast<double>(i);
            const double tick_h = (i == 0 || i == 10) ? 3.0
                                : (i == 5) ? 2.5 : 1.5;  // longer ticks at ends + midpoint for AI orientation
            ModelObject* tick = add_box_object(model,
                                               (boost::format("scale_tick_%1%mm") % i).str(),
                                               tx, scale_bar_y - 0.5 - tick_h / 2.0,
                                               0.45, tick_h, ZCAL_LAYER_THICKNESS);
            set_solid_pad_config(tick);
        }
    }

    // ---- Zone S (solid concentric pad) ----
    {
        ModelObject* mo = add_box_object(model, "zone_S_solid",
                                         zone_S_xc, row_y,
                                         zone, zone, ZCAL_LAYER_THICKNESS);
        set_solid_pad_config(mo);
    }

    // ---- Zone G (3 gap-spacing pairs at 0.5/0.6/0.8mm gaps) ----
    {
        // 3 horizontal pairs of single-perimeter walls. Each pair = 2 walls 20mm long with
        // a specific edge-to-edge gap. Pairs stacked vertically inside the zone, 5mm apart.
        const double wall_len = 20.0;
        const double pair_y_offsets[3] = { +8.0,  0.0, -8.0 };  // top, middle, bottom inside zone
        const double pair_gaps[3]      = {  0.5,  0.6,  0.8 };  // gap to evaluate per pair
        for (int p = 0; p < 3; ++p) {
            const double gap = pair_gaps[p];
            // Two walls centered on (zone_G_xc, row_y + pair_y_offsets[p]); each wall is
            // 0.45mm wide; gap is between their inner edges.
            const double wall_half_gap = (gap + ZCAL_THIN_WALL_WIDTH) / 2.0;
            const double yc_top = row_y + pair_y_offsets[p] + wall_half_gap;
            const double yc_bot = row_y + pair_y_offsets[p] - wall_half_gap;
            ModelObject* w_top = add_box_object(model,
                                                (boost::format("zone_G_gap%1%_top") % gap).str(),
                                                zone_G_xc, yc_top,
                                                wall_len, ZCAL_THIN_WALL_WIDTH, ZCAL_LAYER_THICKNESS);
            ModelObject* w_bot = add_box_object(model,
                                                (boost::format("zone_G_gap%1%_bot") % gap).str(),
                                                zone_G_xc, yc_bot,
                                                wall_len, ZCAL_THIN_WALL_WIDTH, ZCAL_LAYER_THICKNESS);
            set_single_wall_config(w_top);
            set_single_wall_config(w_bot);
        }
    }

    // ---- Zone W (3 freestanding single-perimeter walls) ----
    {
        const double wall_len = 10.0;
        const double wall_dy  = 6.0;  // vertical spacing between walls
        for (int i = 0; i < 3; ++i) {
            const double yc = row_y + (i - 1) * wall_dy;  // -wall_dy, 0, +wall_dy
            ModelObject* mo = add_box_object(model,
                                             (boost::format("zone_W_wall_%1%") % i).str(),
                                             zone_W_xc, yc,
                                             wall_len, ZCAL_THIN_WALL_WIDTH, ZCAL_LAYER_THICKNESS);
            set_single_wall_config(mo);
        }
    }

    // ---- Zone C (4 small corner loops) ----
    {
        struct { double x, y; const char* name; } corners[4] = {
            { -corner_off, -corner_off, "zone_C_BL" },
            { +corner_off, -corner_off, "zone_C_BR" },
            { +corner_off, +corner_off, "zone_C_TR" },
            { -corner_off, +corner_off, "zone_C_TL" },
        };
        const double loop_size = 6.0;  // small ring footprint
        for (auto &c : corners) {
            ModelObject* mo = add_box_object(model, c.name,
                                             c.x, c.y, loop_size, loop_size, ZCAL_LAYER_THICKNESS);
            set_corner_loop_config(mo);
        }
    }

    // ---- Print-config overrides — single layer hard cap ----
    full_config.set_key_value("layer_height", new ConfigOptionFloat(ZCAL_LAYER_THICKNESS));
    full_config.set_key_value("initial_layer_print_height", new ConfigOptionFloat(ZCAL_LAYER_THICKNESS));
    // Disable brim globally (per-object configs handle per-zone behavior; we don't want a
    // global brim wrapping the whole cluster).
    full_config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btNoBrim));
    full_config.set_key_value("brim_width", new ConfigOptionFloat(0.0));
    // Skirt: 1 loop, 4mm offset — useful for the operator to see priming was clean.
    full_config.set_key_value("skirt_loops", new ConfigOptionInt(1));
    full_config.set_key_value("skirt_distance", new ConfigOptionFloat(4.0));
    // Light first-layer speed — the test only matters if the first layer is laid clean.
    full_config.set_key_value("initial_layer_speed", new ConfigOptionFloat(25.0));
    // Resonance/wrapping detection don't apply to a single-layer pad — neutralize both.
    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));
    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_build_zcal_pattern: plate=%1%mm zone=%2%mm fiducials=%3% scale_bar=%4% objects=%5%")
        % plate % zone % params.fiducials % params.scale_bar % model.objects.size();
}

// ----------------------------------------------------------------
// Metadata-comment injection (post-processes the sliced gcode file)
// ----------------------------------------------------------------
//
// Reads the gcode, finds the first non-comment / non-blank line in the header block, and
// prepends a block of coordinate ground-truth comments. AI-vision tooling reads this block
// instead of inferring zone positions from geometry.
bool cli_inject_zcal_metadata(const std::string &gcode_path, const CLIZCalParams &params,
                              const DynamicPrintConfig &full_config)
{
    if (!boost::filesystem::exists(gcode_path)) {
        BOOST_LOG_TRIVIAL(error) << "cli_inject_zcal_metadata: missing gcode file " << gcode_path;
        return false;
    }

    // Resolve bed center from printable_area to convert model coords → bed coords for the metadata.
    Vec2d bed_center(0.0, 0.0);
    if (const ConfigOptionPoints *p = full_config.option<ConfigOptionPoints>("printable_area"); p && !p->values.empty()) {
        BoundingBoxf bb;
        for (const Vec2d &v : p->values) bb.merge(v);
        bed_center = bb.center();
    }

    const double plate    = params.plate_size;
    const double zone     = params.zone_size;
    const double half     = plate / 2.0;
    const double fid_off  = half - 5.0;
    const double corner_off = half - 12.0;
    const double row_y    = 8.0;
    const double zone_pitch = zone + 3.0;
    const double scale_bar_y = -(half - 8.0);

    auto bed_x = [&](double mx) { return bed_center.x() + mx; };
    auto bed_y = [&](double my) { return bed_center.y() + my; };

    std::ostringstream meta;
    meta << "; ---- ORCA Z-OFFSET CALIBRATION METADATA ----\n";
    meta << "; calibration_type = z-offset-pattern\n";
    meta << boost::format("; plate_size_mm = %1$.1f\n") % plate;
    meta << boost::format("; zone_size_mm = %1$.1f\n") % zone;
    meta << boost::format("; bed_center_xy = %1$.2f,%2$.2f\n") % bed_center.x() % bed_center.y();
    if (params.fiducials) {
        meta << boost::format("; fiducial_positions = X%1$.2f,Y%2$.2f X%3$.2f,Y%4$.2f X%5$.2f,Y%6$.2f X%7$.2f,Y%8$.2f\n")
            % bed_x(-fid_off) % bed_y(-fid_off)
            % bed_x(+fid_off) % bed_y(-fid_off)
            % bed_x(+fid_off) % bed_y(+fid_off)
            % bed_x(-fid_off) % bed_y(+fid_off);
    }
    if (params.scale_bar) {
        meta << boost::format("; scale_bar_origin = X%1$.2f,Y%2$.2f\n") % bed_x(-5.0) % bed_y(scale_bar_y);
        meta << "; scale_bar_length_mm = 10\n";
        meta << "; scale_bar_tick_count = 11\n";
    }
    meta << boost::format("; zone_S_center = X%1$.2f,Y%2$.2f  zone_S_size = %3$.1f\n") % bed_x(-zone_pitch) % bed_y(row_y) % zone;
    meta << boost::format("; zone_G_center = X%1$.2f,Y%2$.2f  zone_G_size = %3$.1f  zone_G_gaps_mm = 0.5,0.6,0.8\n") % bed_x(0.0) % bed_y(row_y) % zone;
    meta << boost::format("; zone_W_center = X%1$.2f,Y%2$.2f  zone_W_size = %3$.1f  zone_W_wall_count = 3\n") % bed_x(+zone_pitch) % bed_y(row_y) % zone;
    meta << boost::format("; zone_C_corners = X%1$.2f,Y%2$.2f X%3$.2f,Y%4$.2f X%5$.2f,Y%6$.2f X%7$.2f,Y%8$.2f\n")
        % bed_x(-corner_off) % bed_y(-corner_off)
        % bed_x(+corner_off) % bed_y(-corner_off)
        % bed_x(+corner_off) % bed_y(+corner_off)
        % bed_x(-corner_off) % bed_y(+corner_off);
    meta << "; ---- END Z-OFFSET CALIBRATION METADATA ----\n";

    // Insert right after the HEADER_BLOCK_START marker (consistent with other Orca metadata).
    std::ifstream in(gcode_path);
    if (!in) return false;
    std::ostringstream body;
    body << in.rdbuf();
    in.close();
    std::string contents = body.str();

    const std::string marker = "; HEADER_BLOCK_START";
    const size_t marker_pos = contents.find(marker);
    size_t insert_pos = 0;
    if (marker_pos != std::string::npos) {
        const size_t line_end = contents.find('\n', marker_pos);
        insert_pos = (line_end != std::string::npos) ? line_end + 1 : marker_pos + marker.size();
    } else {
        // No HEADER_BLOCK marker — fall through to the very top, after any leading shebang.
        insert_pos = 0;
    }
    contents.insert(insert_pos, meta.str());

    std::ofstream out(gcode_path, std::ios::trunc);
    if (!out) return false;
    out.write(contents.data(), contents.size());
    out.close();
    BOOST_LOG_TRIVIAL(info) << "cli_inject_zcal_metadata: wrote metadata block to " << gcode_path;
    return true;
}

}
