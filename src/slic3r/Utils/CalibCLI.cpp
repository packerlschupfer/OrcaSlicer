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

#include "libslic3r/CutUtils.hpp"
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
    if (name == "temp-tower")               return CLICalibType::TempTower;
    if (name == "vol-speed-tower")          return CLICalibType::VolSpeedTower;
    if (name == "pa-tower")                 return CLICalibType::PATower;
    if (name == "retraction-tower")         return CLICalibType::RetractionTower;
    if (name == "vfa-tower")                return CLICalibType::VFATower;
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
    case CLICalibType::TempTower:
        return root + "/calib/temperature_tower/temperature_tower.drc";
    case CLICalibType::VolSpeedTower:
        return root + "/calib/volumetric_speed/SpeedTestStructure.drc";
    case CLICalibType::PATower:
        return root + "/calib/pressure_advance/tower_with_seam.drc";
    case CLICalibType::RetractionTower:
        return root + "/calib/retraction/retraction_tower.drc";
    case CLICalibType::VFATower:
        return root + "/calib/vfa/vfa.drc";
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

    // ---- Peelable structural frame (4 sides of a ring around the fiducials) ----
    //      Multi-layer so the whole plate peels off the textured PEI sheet as one piece without
    //      tearing the single-layer calibration surface. Frame is OUTSIDE the fiducial bbox by
    //      params.frame_margin, frame_width thick, frame_layers tall. Each side is a solid
    //      filled box; per-object config gives it strong walls + dense top so it survives peeling.
    if (params.frame && params.frame_layers > 0 && params.frame_width > 0.0) {
        // fiducial outer edge in model space: |fid_off| + (fid_size/2)
        const double fid_outer        = fid_off + (5.0 / 2.0);
        const double frame_inner_edge = fid_outer + params.frame_margin;
        const double frame_outer_edge = frame_inner_edge + params.frame_width;
        const double frame_height     = ZCAL_LAYER_THICKNESS * static_cast<double>(params.frame_layers);
        // Horizontal bars span the full outer width; vertical bars sit between them to avoid
        // double-extrusion in the corners.
        const double horiz_len = 2.0 * frame_outer_edge;
        const double vert_len  = 2.0 * frame_inner_edge;
        const double bar_xc    = (frame_inner_edge + frame_outer_edge) / 2.0;
        const double bar_yc    = (frame_inner_edge + frame_outer_edge) / 2.0;

        ModelObject* top    = add_box_object(model, "frame_top",    0.0,        +bar_yc, horiz_len, params.frame_width, frame_height);
        ModelObject* bottom = add_box_object(model, "frame_bottom", 0.0,        -bar_yc, horiz_len, params.frame_width, frame_height);
        ModelObject* left   = add_box_object(model, "frame_left",   -bar_xc,    0.0,     params.frame_width, vert_len, frame_height);
        ModelObject* right  = add_box_object(model, "frame_right",  +bar_xc,    0.0,     params.frame_width, vert_len, frame_height);

        // Strong walls + dense top so the multi-layer bar survives being peeled. With a 3mm-wide
        // box and 0.45mm line width, 4 walls = 1.8mm and the remainder fills with concentric.
        for (ModelObject* fb : {top, bottom, left, right}) {
            fb->config.set_key_value("wall_loops", new ConfigOptionInt(3));
            fb->config.set_key_value("top_shell_layers", new ConfigOptionInt(std::min(params.frame_layers, 2)));
            fb->config.set_key_value("bottom_shell_layers", new ConfigOptionInt(std::min(params.frame_layers, 2)));
            fb->config.set_key_value("sparse_infill_density", new ConfigOptionPercent(100));
            fb->config.set_key_value("top_surface_pattern", new ConfigOptionEnum<InfillPattern>(ipConcentric));
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


// ============================================================
// Tower-shaped calibration tests (Temperature, Volumetric Speed, PA Tower).
// ============================================================

namespace {

// Headless port of Plater::cut_horizontal (Plater.cpp:12877). The GUI uses Cut + apply_cut_object_to_model
// which clears canvas + obj-list state we don't have on the CLI side. Here: build the cut, replace the
// original object in model.objects with the cut result (KeepLower drops the top, KeepUpper drops the
// bottom — matching the GUI's KeepLower/KeepUpper semantics).
void cli_cut_horizontal(Model &model, size_t obj_idx, double z, ModelObjectCutAttributes attributes)
{
    if (obj_idx >= model.objects.size()) {
        BOOST_LOG_TRIVIAL(error) << "cli_cut_horizontal: obj_idx " << obj_idx << " out of range";
        return;
    }
    if (!attributes.has(ModelObjectCutAttribute::KeepUpper) && !attributes.has(ModelObjectCutAttribute::KeepLower))
        return;

    ModelObject *src = model.objects[obj_idx];
    if (src->instances.empty()) return;

    const Vec3d instance_offset = src->instances[0]->get_offset();
    Cut cut(src, /*instance_idx=*/0,
            Geometry::translation_transform(z * Vec3d::UnitZ() - instance_offset),
            attributes);
    const auto &new_objects = cut.perform_with_plane();
    if (new_objects.empty()) {
        BOOST_LOG_TRIVIAL(warning) << "cli_cut_horizontal: perform_with_plane produced no objects";
        return;
    }

    // Replace the original object with the cut result(s). KeepLower-only typically yields a single
    // ModelObject; KeepUpper-only same. KeepBoth would yield two.
    model.delete_object(obj_idx);
    for (ModelObject *new_obj : new_objects)
        model.add_object(*new_obj);
}

// First-extruder nozzle diameter helper. All three towers use slot 0.
double cli_nozzle_diameter(const DynamicPrintConfig &full_config, double fallback = 0.4)
{
    const ConfigOptionFloats *nd = full_config.option<ConfigOptionFloats>("nozzle_diameter");
    if (nd && !nd->values.empty()) return nd->values[0];
    return fallback;
}

} // anonymous

//ORCA: Headless port of Plater::calib_temp (Plater.cpp:13178). Loads the temperature_tower mesh,
//      cuts off the unused top + bottom based on params.start/end (with the 5°C step the GUI
//      hard-codes), scales by nozzle ratio, and applies the per-object + filament + print configs
//      line-for-line. Calls below need set_calib_params() on the Print engine before slicing.
void cli_apply_temp_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params)
{
    if (model.objects.empty()) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_temp_tower: model is empty";
        return;
    }
    constexpr double base_temp_tower_nozzle_diameter = 0.4;
    constexpr double base_temp_tower_block_height    = 10.0;
    constexpr int    base_temp_tower_temp_step       = 5;

    const long start_temp = std::lround(params.start);

    const double nozzle_diameter = std::max(cli_nozzle_diameter(full_config, base_temp_tower_nozzle_diameter), 0.001);
    const double nozzle_scale    = nozzle_diameter / base_temp_tower_nozzle_diameter;

    // Cut upper (drop blocks colder than params.end).
    {
        const auto obj_bb = model.objects[0]->bounding_box_exact();
        const long block_count = std::lround((500 - params.end) / base_temp_tower_temp_step + 1);
        if (block_count > 0) {
            const double new_height = block_count * base_temp_tower_block_height - EPSILON;
            if (new_height < obj_bb.size().z())
                cli_cut_horizontal(model, 0, new_height, ModelObjectCutAttribute::KeepLower);
        }
    }
    // Cut bottom (drop blocks hotter than params.start).
    {
        const auto obj_bb = model.objects[0]->bounding_box_exact();
        const long block_count = std::lround((500 - params.start) / base_temp_tower_temp_step);
        if (block_count > 0) {
            const double new_height = block_count * base_temp_tower_block_height + EPSILON;
            if (new_height < obj_bb.size().z())
                cli_cut_horizontal(model, 0, new_height, ModelObjectCutAttribute::KeepUpper);
        }
    }

    if (std::abs(nozzle_scale - 1.0) > EPSILON)
        model.objects[0]->scale(nozzle_scale, nozzle_scale, nozzle_scale);
    model.objects[0]->ensure_on_bed();

    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));
    full_config.set_key_value("nozzle_temperature_initial_layer", new ConfigOptionInts(1, static_cast<int>(start_temp)));
    full_config.set_key_value("nozzle_temperature", new ConfigOptionInts(1, static_cast<int>(start_temp)));
    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));
    full_config.set_key_value("initial_layer_print_height", new ConfigOptionFloat(nozzle_diameter / 2));

    auto &obj_cfg = model.objects[0]->config;
    obj_cfg.set_key_value("layer_height", new ConfigOptionFloat(nozzle_diameter / 2));
    obj_cfg.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btOuterOnly));
    obj_cfg.set_key_value("brim_width", new ConfigOptionFloat(5.0));
    obj_cfg.set_key_value("brim_object_gap", new ConfigOptionFloat(0.0));
    obj_cfg.set_key_value("alternate_extra_wall", new ConfigOptionBool(false));
    obj_cfg.set_key_value("seam_slope_type", new ConfigOptionEnum<SeamScarfType>(SeamScarfType::None));
    obj_cfg.set_key_value("overhang_reverse", new ConfigOptionBool(false));
    obj_cfg.set_key_value("precise_z_height", new ConfigOptionBool(false));

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_temp_tower: start=%1%°C end=%2%°C nozzle=%3$.2fmm scale=%4$.3f")
        % params.start % params.end % nozzle_diameter % nozzle_scale;
}

//ORCA: Headless port of Plater::calib_max_vol_speed (Plater.cpp:13258). Loads SpeedTestStructure,
//      scales horizontally to fit the bed, cuts the upper to match the user's flow range, applies
//      configs line-for-line. The Print engine's calib_params is set to speed (mm/s), not flow (mm³/s):
//      that conversion is done by cli_tower_get_calib_params() before set_calib_params().
void cli_apply_vol_speed_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params)
{
    if (model.objects.empty()) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_vol_speed_tower: model is empty";
        return;
    }
    auto *obj = model.objects[0];

    // Horizontal scale so the tower fits on the bed minus 10mm margin.
    if (const auto *area = full_config.option<ConfigOptionPoints>("printable_area"); area && area->values.size() >= 4) {
        BoundingBoxf bed_ext = get_extents(area->values);
        const double scale_obj = (bed_ext.size().x() - 10.0) / obj->bounding_box_exact().size().x();
        if (scale_obj < 1.0)
            obj->scale(scale_obj, 1.0, 1.0);
    }

    const double nozzle_diameter = cli_nozzle_diameter(full_config);
    const double line_width      = nozzle_diameter * 1.75;
    const double layer_height    = nozzle_diameter * 0.8;

    if (auto *max_lh = full_config.option<ConfigOptionFloats>("max_layer_height"); max_lh && !max_lh->values.empty()) {
        if (max_lh->values[0] < layer_height) max_lh->values[0] = layer_height;
    }

    full_config.set_key_value("filament_max_volumetric_speed", new ConfigOptionFloats{200.0});
    full_config.set_key_value("slow_down_layer_time", new ConfigOptionFloats{0.0});
    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));
    full_config.set_key_value("timelapse_type", new ConfigOptionEnum<TimelapseType>(tlTraditional));
    full_config.set_key_value("spiral_mode", new ConfigOptionBool(true));
    full_config.set_key_value("max_volumetric_extrusion_rate_slope", new ConfigOptionFloat(0));
    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));

    auto &obj_cfg = obj->config;
    obj_cfg.set_key_value("enable_overhang_speed", new ConfigOptionBool(false));
    obj_cfg.set_key_value("wall_loops", new ConfigOptionInt(1));
    obj_cfg.set_key_value("alternate_extra_wall", new ConfigOptionBool(false));
    obj_cfg.set_key_value("top_shell_layers", new ConfigOptionInt(0));
    obj_cfg.set_key_value("bottom_shell_layers", new ConfigOptionInt(0));
    obj_cfg.set_key_value("sparse_infill_density", new ConfigOptionPercent(0));
    obj_cfg.set_key_value("outer_wall_line_width", new ConfigOptionFloatOrPercent(line_width, false));
    obj_cfg.set_key_value("layer_height", new ConfigOptionFloat(layer_height));
    obj_cfg.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btOuterAndInner));
    obj_cfg.set_key_value("brim_width", new ConfigOptionFloat(5.0));
    obj_cfg.set_key_value("brim_object_gap", new ConfigOptionFloat(0.0));
    obj_cfg.set_key_value("precise_z_height", new ConfigOptionBool(false));

    // Cut upper to the height that covers the user's range.
    if (params.step > 0.0) {
        const double obj_z   = obj->bounding_box_exact().size().z();
        const double height  = (params.end - params.start + 1.0) / params.step;
        if (height < obj_z)
            cli_cut_horizontal(model, 0, height, ModelObjectCutAttribute::KeepLower);
    }

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_vol_speed_tower: start=%1% end=%2% step=%3% mm³/s, line_width=%4$.3f layer_h=%5$.3f")
        % params.start % params.end % params.step % line_width % layer_height;
}

//ORCA: Headless port of Plater::_calib_pa_tower (Plater.cpp:12896). Loads the PA tower with seam,
//      cuts upper based on (end-start)/step, applies the configs line-for-line. The find_optimal_PA_speed
//      call uses CalibPressureAdvance::find_optimal_PA_speed — same helper the GUI uses.
void cli_apply_pa_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params)
{
    if (model.objects.empty()) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_pa_tower: model is empty";
        return;
    }
    const double nozzle_diameter = cli_nozzle_diameter(full_config);

    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));
    full_config.set_key_value("slow_down_layer_time", new ConfigOptionFloats{1.0});
    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));
    full_config.set_key_value("overhang_reverse", new ConfigOptionBool(false));
    full_config.set_key_value("precise_z_height", new ConfigOptionBool(false));
    full_config.set_key_value("max_volumetric_extrusion_rate_slope", new ConfigOptionFloat(0));

    auto &obj_cfg = model.objects[0]->config;
    obj_cfg.set_key_value("alternate_extra_wall", new ConfigOptionBool(false));

    // GUI uses CalibPressureAdvance::find_optimal_PA_speed — reuse it verbatim.
    const double line_width   = full_config.get_abs_value("line_width", nozzle_diameter);
    const double layer_height = full_config.get_abs_value("layer_height");
    const double wall_speed   = CalibPressureAdvance::find_optimal_PA_speed(full_config, line_width, layer_height, 0, 0);
    obj_cfg.set_key_value("outer_wall_speed", new ConfigOptionFloat(wall_speed));
    obj_cfg.set_key_value("inner_wall_speed", new ConfigOptionFloat(wall_speed));
    obj_cfg.set_key_value("seam_position", new ConfigOptionEnum<SeamPosition>(spRear));
    obj_cfg.set_key_value("wall_loops", new ConfigOptionInt(2));
    obj_cfg.set_key_value("top_shell_layers", new ConfigOptionInt(0));
    obj_cfg.set_key_value("bottom_shell_layers", new ConfigOptionInt(0));
    obj_cfg.set_key_value("sparse_infill_density", new ConfigOptionPercent(0));
    obj_cfg.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btEar));
    obj_cfg.set_key_value("brim_object_gap", new ConfigOptionFloat(0.0));
    obj_cfg.set_key_value("brim_ears_max_angle", new ConfigOptionFloat(135.0));
    obj_cfg.set_key_value("brim_width", new ConfigOptionFloat(6.0));
    obj_cfg.set_key_value("seam_slope_type", new ConfigOptionEnum<SeamScarfType>(SeamScarfType::None));

    if (params.step > 0.0) {
        const double new_height = std::ceil((params.end - params.start) / params.step) + 1.0;
        const auto obj_bb = model.objects[0]->bounding_box_exact();
        if (new_height < obj_bb.size().z())
            cli_cut_horizontal(model, 0, new_height, ModelObjectCutAttribute::KeepLower);
    }

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_pa_tower: start=%1% end=%2% step=%3% s, wall_speed=%4$.2fmm/s")
        % params.start % params.end % params.step % wall_speed;
}

//ORCA: Headless port of Plater::calib_retraction (Plater.cpp:13333). Loads retraction_tower.drc,
//      cuts to user's range, applies per-object configs + a layer-height pick driven by nozzle
//      size. The Print engine's per-layer modulation handles the retraction-length progression.
void cli_apply_retraction_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params)
{
    if (model.objects.empty()) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_retraction_tower: model is empty";
        return;
    }
    auto *obj = model.objects[0];

    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));

    const double nozzle_diameter = cli_nozzle_diameter(full_config);
    double layer_height = 0.2;
    if      (nozzle_diameter <= 0.1) layer_height = 0.05;
    else if (nozzle_diameter <= 0.2) layer_height = 0.10;

    if (auto *max_lh = full_config.option<ConfigOptionFloats>("max_layer_height"); max_lh && !max_lh->values.empty())
        if (max_lh->values[0] < layer_height) max_lh->values[0] = layer_height;

    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));
    full_config.set_key_value("use_firmware_retraction", new ConfigOptionBool(false));
    full_config.set_key_value("initial_layer_print_height", new ConfigOptionFloat(layer_height));

    auto &cfg = obj->config;
    cfg.set_key_value("wall_loops", new ConfigOptionInt(2));
    cfg.set_key_value("top_shell_layers", new ConfigOptionInt(0));
    cfg.set_key_value("bottom_shell_layers", new ConfigOptionInt(3));
    cfg.set_key_value("sparse_infill_density", new ConfigOptionPercent(0));
    cfg.set_key_value("layer_height", new ConfigOptionFloat(layer_height));
    cfg.set_key_value("alternate_extra_wall", new ConfigOptionBool(false));
    cfg.set_key_value("seam_position", new ConfigOptionEnum<SeamPosition>(spAligned));
    cfg.set_key_value("wall_sequence", new ConfigOptionEnum<WallSequence>(WallSequence::InnerOuter));
    cfg.set_key_value("overhang_reverse", new ConfigOptionBool(false));
    cfg.set_key_value("precise_z_height", new ConfigOptionBool(false));

    if (params.step > 0.0) {
        const auto obj_bb = obj->bounding_box_exact();
        const double height = 1.0 + 0.4 + (params.end - params.start) / params.step - EPSILON;
        if (height < obj_bb.size().z())
            cli_cut_horizontal(model, 0, height, ModelObjectCutAttribute::KeepLower);
    }

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_retraction_tower: start=%1% end=%2% step=%3% mm, layer_h=%4$.3f")
        % params.start % params.end % params.step % layer_height;
}

//ORCA: Headless port of Plater::calib_VFA (Plater.cpp:13391). Loads vfa.drc, applies spiral-mode
//      single-perimeter configs, cuts to range. The Print engine modulates outer_wall_speed at
//      5mm Z intervals (GCode.cpp:4612).
void cli_apply_vfa_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params)
{
    if (model.objects.empty()) {
        BOOST_LOG_TRIVIAL(error) << "cli_apply_vfa_tower: model is empty";
        return;
    }
    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));
    full_config.set_key_value("slow_down_layer_time", new ConfigOptionFloats{0.0});
    full_config.set_key_value("enable_overhang_speed", new ConfigOptionBool(false));
    full_config.set_key_value("timelapse_type", new ConfigOptionEnum<TimelapseType>(tlTraditional));
    full_config.set_key_value("wall_loops", new ConfigOptionInt(1));
    full_config.set_key_value("alternate_extra_wall", new ConfigOptionBool(false));
    full_config.set_key_value("top_shell_layers", new ConfigOptionInt(0));
    full_config.set_key_value("bottom_shell_layers", new ConfigOptionInt(1));
    full_config.set_key_value("sparse_infill_density", new ConfigOptionPercent(0));
    full_config.set_key_value("detect_thin_wall", new ConfigOptionBool(false));
    full_config.set_key_value("spiral_mode", new ConfigOptionBool(true));
    full_config.set_key_value("enable_wrapping_detection", new ConfigOptionBool(false));
    full_config.set_key_value("precise_z_height", new ConfigOptionBool(false));

    auto &cfg = model.objects[0]->config;
    cfg.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btOuterOnly));
    cfg.set_key_value("brim_width", new ConfigOptionFloat(3.0));
    cfg.set_key_value("brim_object_gap", new ConfigOptionFloat(0.0));

    if (params.step > 0.0) {
        const auto obj_bb = model.objects[0]->bounding_box_exact();
        const double height = 5.0 * ((params.end - params.start) / params.step + 1.0);
        if (height < obj_bb.size().z())
            cli_cut_horizontal(model, 0, height, ModelObjectCutAttribute::KeepLower);
    }

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_vfa_tower: start=%1% end=%2% step=%3% mm/s")
        % params.start % params.end % params.step;
}

bool cli_tower_get_calib_params(CLICalibType type, const CLITowerParams &params,
                                const DynamicPrintConfig &full_config, Calib_Params &out)
{
    out = Calib_Params{};
    out.start = params.start;
    out.end   = params.end;
    out.step  = params.step;

    switch (type) {
    case CLICalibType::TempTower:
        out.mode = CalibMode::Calib_Temp_Tower;
        return true;
    case CLICalibType::VolSpeedTower: {
        out.mode = CalibMode::Calib_Vol_speed_Tower;
        // Convert user mm³/s to print-engine mm/s using mm³_per_mm = Flow * filament_flow_ratio (slot 0).
        const double nozzle_diameter = cli_nozzle_diameter(full_config);
        const double line_width      = nozzle_diameter * 1.75;
        const double layer_height    = nozzle_diameter * 0.8;
        double flow_ratio = 1.0;
        if (auto *fr = full_config.option<ConfigOptionFloatsNullable>("filament_flow_ratio"); fr && !fr->values.empty())
            flow_ratio = fr->values[0];
        const double mm3_per_mm = Flow(line_width, layer_height, nozzle_diameter).mm3_per_mm() * flow_ratio;
        if (mm3_per_mm > 0.0) {
            out.start = params.start / mm3_per_mm;
            out.end   = params.end   / mm3_per_mm;
            out.step  = params.step  / mm3_per_mm;
        }
        return true;
    }
    case CLICalibType::PATower:
        out.mode = CalibMode::Calib_PA_Tower;
        return true;
    case CLICalibType::RetractionTower:
        out.mode = CalibMode::Calib_Retraction_tower;
        return true;
    case CLICalibType::VFATower:
        out.mode = CalibMode::Calib_VFA_Tower;
        return true;
    default:
        return false;
    }
}

// ============================================================
// Unified calibration outputs (gcode header metadata + JSON sidecar).
// ============================================================

namespace {

// Tiny JSON writer — std::string concat. The structure is fixed enough that pulling in nlohmann/json
// or boost::json would be overkill, and we already control every key/value being emitted.
struct JsonOut {
    std::string body;
    int indent = 0;
    void w(const std::string &s) { body += std::string(indent * 2, ' '); body += s; }
    static std::string str(const std::string &s) {
        std::string r = "\"";
        for (char c : s) {
            if (c == '"' || c == '\\') { r += '\\'; r += c; }
            else if (c == '\n') r += "\\n";
            else r += c;
        }
        r += '"';
        return r;
    }
    static std::string num(double d) {
        std::ostringstream os;
        os.precision(6);
        os << std::fixed << d;
        std::string s = os.str();
        // Trim trailing zeros after the decimal point for readability.
        if (s.find('.') != std::string::npos) {
            while (!s.empty() && s.back() == '0') s.pop_back();
            if (!s.empty() && s.back() == '.') s.pop_back();
        }
        return s;
    }
};

Vec2d resolve_bed_center(const DynamicPrintConfig &full_config)
{
    if (const ConfigOptionPoints *p = full_config.option<ConfigOptionPoints>("printable_area"); p && !p->values.empty()) {
        BoundingBoxf bb;
        for (const Vec2d &v : p->values) bb.merge(v);
        return bb.center();
    }
    return Vec2d{125.0, 110.0};
}

// Returns the calibration type name as a stable string (mirrors the CLI --calibrate-type values).
const char* calib_type_name(CLICalibType t)
{
    switch (t) {
    case CLICalibType::FlowRate_YOLO_Recommended:    return "flow-yolo-recommended";
    case CLICalibType::FlowRate_YOLO_Perfectionist:  return "flow-yolo-perfectionist";
    case CLICalibType::FlowRate_Pass1:               return "flow-pass1";
    case CLICalibType::FlowRate_Pass2:               return "flow-pass2";
    case CLICalibType::ZOffsetPattern:               return "z-offset-pattern";
    case CLICalibType::TempTower:                    return "temp-tower";
    case CLICalibType::VolSpeedTower:                return "vol-speed-tower";
    case CLICalibType::PATower:                      return "pa-tower";
    case CLICalibType::RetractionTower:              return "retraction-tower";
    case CLICalibType::VFATower:                     return "vfa-tower";
    default:                                         return "none";
    }
}

// For tower types, emit a per-block table of {z_low, z_high, value} so AI can look up
// "scan defect at Z=27mm" → block params without inverting the start/end/step math.
// Encodes the Print engine's per-layer modulation formula from GCode.cpp:4604-4630.
struct TowerBlock { double z_low, z_high, value; };
std::vector<TowerBlock> compute_tower_blocks(CLICalibType type, const CLITowerParams &p)
{
    std::vector<TowerBlock> out;
    if (p.step <= 0.0) return out;

    auto push = [&](double zl, double zh, double v) { out.push_back({zl, zh, v}); };

    switch (type) {
    case CLICalibType::TempTower: {
        // 10mm per block, 5°C step, start at bottom (hottest) → end at top (coolest).
        const int blocks = static_cast<int>(std::round((p.start - p.end) / 5.0) + 1);
        for (int i = 0; i < blocks; ++i)
            push(i * 10.0, (i + 1) * 10.0, p.start - i * 5.0);
        break;
    }
    case CLICalibType::PATower: {
        // GCode.cpp:4605: PA = start + (int)Z * step. So blocks of 1mm height.
        const int blocks = static_cast<int>(std::round((p.end - p.start) / p.step) + 1);
        for (int i = 0; i < blocks; ++i)
            push(static_cast<double>(i), i + 1.0, p.start + i * p.step);
        break;
    }
    case CLICalibType::VolSpeedTower: {
        // GCode.cpp:4618: speed_mm_per_s = start + Z * step (after mm3/s→mm/s conversion).
        // For user-facing metadata we emit the ORIGINAL mm3/s values, blockless (Z-range = full cut).
        const int steps = static_cast<int>(std::round((p.end - p.start) / p.step) + 1);
        // height-per-step depends on the conversion; we just expose the linear formula.
        for (int i = 0; i < steps; ++i)
            push(static_cast<double>(i), i + 1.0, p.start + i * p.step);
        break;
    }
    case CLICalibType::RetractionTower: {
        // GCode.cpp:4623: length = start + floor(max(0,Z-0.4)) * step. 1mm-per-block.
        const int blocks = static_cast<int>(std::round((p.end - p.start) / p.step) + 1);
        for (int i = 0; i < blocks; ++i)
            push(0.4 + i, 0.4 + i + 1.0, p.start + i * p.step);
        break;
    }
    case CLICalibType::VFATower: {
        // GCode.cpp:4613: speed = start + floor(Z / 5.0) * step. 5mm-per-block.
        const int blocks = static_cast<int>(std::round((p.end - p.start) / p.step) + 1);
        for (int i = 0; i < blocks; ++i)
            push(i * 5.0, (i + 1) * 5.0, p.start + i * p.step);
        break;
    }
    default: break;
    }
    return out;
}

const char* tower_value_key(CLICalibType t)
{
    switch (t) {
    case CLICalibType::TempTower:        return "temp_c";
    case CLICalibType::PATower:          return "pa_s";
    case CLICalibType::VolSpeedTower:    return "vol_speed_mm3_per_s";
    case CLICalibType::RetractionTower:  return "retraction_mm";
    case CLICalibType::VFATower:         return "speed_mm_per_s";
    default: return "value";
    }
}

// Compute per-block info for flow-rate tests by walking model.objects and parsing names.
struct FlowBlockInfo { std::string name; double modifier; double print_flow_ratio; double x, y; };
std::vector<FlowBlockInfo> compute_flow_blocks(const Model *model, double cur_flowrate, bool is_linear)
{
    std::vector<FlowBlockInfo> out;
    if (!model) return out;
    for (const ModelObject *mo : model->objects) {
        if (!mo) continue;
        bool ok = false;
        const float mod = parse_flow_modifier(mo->name, ok);
        if (!ok) continue;
        const double pfr = is_linear
            ? (cur_flowrate + mod) / cur_flowrate
            : 1.0 + mod / 100.0;
        Vec3d off = mo->instances.empty() ? Vec3d::Zero() : mo->instances[0]->get_offset();
        out.push_back({mo->name, mod, pfr, off.x(), off.y()});
    }
    return out;
}

} // anonymous

bool cli_emit_calib_outputs(const std::string &gcode_path,
                            CLICalibType type,
                            const DynamicPrintConfig &full_config,
                            const CLIZCalParams *zcal_params,
                            const CLITowerParams *tower_params,
                            const CLIFlowRateParams *flow_params,
                            const Model *model)
{
    if (type == CLICalibType::NoCalib) return false;
    if (!boost::filesystem::exists(gcode_path)) {
        BOOST_LOG_TRIVIAL(error) << "cli_emit_calib_outputs: missing gcode file " << gcode_path;
        return false;
    }

    const Vec2d bed_center = resolve_bed_center(full_config);
    const std::string type_name = calib_type_name(type);

    // ---- Compose gcode header metadata block (machine-readable comments) ----
    std::ostringstream meta;
    meta << "; ---- ORCA CALIBRATION METADATA ----\n";
    meta << "; calibration_type = " << type_name << "\n";
    meta << boost::format("; bed_center_xy = %1$.2f,%2$.2f\n") % bed_center.x() % bed_center.y();

    // ---- Compose JSON sidecar payload ----
    JsonOut J;
    J.indent = 0;
    J.w("{\n"); J.indent = 1;
    J.w(JsonOut::str("calibration_type") + ": " + JsonOut::str(type_name) + ",\n");
    J.w(JsonOut::str("bed_center") + ": {"
        + JsonOut::str("x") + ": " + JsonOut::num(bed_center.x()) + ", "
        + JsonOut::str("y") + ": " + JsonOut::num(bed_center.y()) + "},\n");
    J.w(JsonOut::str("gcode_file") + ": " + JsonOut::str(boost::filesystem::path(gcode_path).filename().string()) + ",\n");

    if (type == CLICalibType::ZOffsetPattern && zcal_params) {
        const auto &p = *zcal_params;
        const double half       = p.plate_size / 2.0;
        const double fid_off    = half - 5.0;
        const double corner_off = half - 12.0;
        const double row_y      = 8.0;
        const double zone_pitch = p.zone_size + 3.0;
        const double scale_bar_y = -(half - 8.0);

        meta << boost::format("; plate_size_mm = %1$.1f\n") % p.plate_size;
        meta << boost::format("; zone_size_mm = %1$.1f\n") % p.zone_size;
        if (p.fiducials) {
            meta << boost::format("; fiducial_positions = X%1$.2f,Y%2$.2f X%3$.2f,Y%4$.2f X%5$.2f,Y%6$.2f X%7$.2f,Y%8$.2f\n")
                % (bed_center.x() - fid_off) % (bed_center.y() - fid_off)
                % (bed_center.x() + fid_off) % (bed_center.y() - fid_off)
                % (bed_center.x() + fid_off) % (bed_center.y() + fid_off)
                % (bed_center.x() - fid_off) % (bed_center.y() + fid_off);
        }
        if (p.scale_bar) {
            meta << boost::format("; scale_bar_origin = X%1$.2f,Y%2$.2f\n") % (bed_center.x() - 5.0) % (bed_center.y() + scale_bar_y);
            meta << "; scale_bar_length_mm = 10\n";
            meta << "; scale_bar_tick_count = 11\n";
        }
        meta << boost::format("; zone_S_center = X%1$.2f,Y%2$.2f  zone_S_size = %3$.1f\n")
            % (bed_center.x() - zone_pitch) % (bed_center.y() + row_y) % p.zone_size;
        meta << boost::format("; zone_G_center = X%1$.2f,Y%2$.2f  zone_G_size = %3$.1f  zone_G_gaps_mm = 0.5,0.6,0.8\n")
            % bed_center.x() % (bed_center.y() + row_y) % p.zone_size;
        meta << boost::format("; zone_W_center = X%1$.2f,Y%2$.2f  zone_W_size = %3$.1f  zone_W_wall_count = 3\n")
            % (bed_center.x() + zone_pitch) % (bed_center.y() + row_y) % p.zone_size;
        meta << boost::format("; zone_C_corners = X%1$.2f,Y%2$.2f X%3$.2f,Y%4$.2f X%5$.2f,Y%6$.2f X%7$.2f,Y%8$.2f\n")
            % (bed_center.x() - corner_off) % (bed_center.y() - corner_off)
            % (bed_center.x() + corner_off) % (bed_center.y() - corner_off)
            % (bed_center.x() + corner_off) % (bed_center.y() + corner_off)
            % (bed_center.x() - corner_off) % (bed_center.y() + corner_off);

        // JSON
        J.w(JsonOut::str("plate_size_mm") + ": " + JsonOut::num(p.plate_size) + ",\n");
        J.w(JsonOut::str("zone_size_mm") + ": " + JsonOut::num(p.zone_size) + ",\n");
        J.w(JsonOut::str("zones") + ": [\n"); J.indent = 2;
        auto zone = [&](const std::string &name, double mx, double my, double sz, const std::string &extra) {
            J.w("{" + JsonOut::str("name") + ": " + JsonOut::str(name)
                + ", " + JsonOut::str("center") + ": {"
                + JsonOut::str("x") + ": " + JsonOut::num(bed_center.x() + mx) + ", "
                + JsonOut::str("y") + ": " + JsonOut::num(bed_center.y() + my) + "}"
                + ", " + JsonOut::str("size_mm") + ": " + JsonOut::num(sz)
                + (extra.empty() ? std::string{} : (", " + extra))
                + "}");
        };
        zone("S", -zone_pitch, row_y, p.zone_size, JsonOut::str("kind") + ": " + JsonOut::str("solid_concentric"));
        J.body += ",\n";
        zone("G",  0.0,        row_y, p.zone_size, JsonOut::str("kind") + ": " + JsonOut::str("gap_pairs") + ", " + JsonOut::str("gaps_mm") + ": [0.5,0.6,0.8]");
        J.body += ",\n";
        zone("W", +zone_pitch, row_y, p.zone_size, JsonOut::str("kind") + ": " + JsonOut::str("single_walls") + ", " + JsonOut::str("count") + ": 3");
        J.body += "\n"; J.indent = 1;
        J.w("],\n");
        if (p.fiducials) {
            J.w(JsonOut::str("fiducials") + ": [");
            const double xs[4] = {bed_center.x() - fid_off, bed_center.x() + fid_off, bed_center.x() + fid_off, bed_center.x() - fid_off};
            const double ys[4] = {bed_center.y() - fid_off, bed_center.y() - fid_off, bed_center.y() + fid_off, bed_center.y() + fid_off};
            for (int i = 0; i < 4; ++i) {
                if (i) J.body += ",";
                J.body += "{" + JsonOut::str("x") + ":" + JsonOut::num(xs[i]) + "," + JsonOut::str("y") + ":" + JsonOut::num(ys[i]) + "}";
            }
            J.body += "],\n";
        }
        if (p.scale_bar) {
            J.w(JsonOut::str("scale_bar") + ": {"
                + JsonOut::str("origin") + ":{" + JsonOut::str("x") + ":" + JsonOut::num(bed_center.x() - 5.0)
                + "," + JsonOut::str("y") + ":" + JsonOut::num(bed_center.y() + scale_bar_y) + "},"
                + JsonOut::str("length_mm") + ":10,"
                + JsonOut::str("tick_count") + ":11},\n");
        }
        J.w(JsonOut::str("corner_loops_C") + ": [");
        const double cxs[4] = {bed_center.x() - corner_off, bed_center.x() + corner_off, bed_center.x() + corner_off, bed_center.x() - corner_off};
        const double cys[4] = {bed_center.y() - corner_off, bed_center.y() - corner_off, bed_center.y() + corner_off, bed_center.y() + corner_off};
        for (int i = 0; i < 4; ++i) {
            if (i) J.body += ",";
            J.body += "{" + JsonOut::str("x") + ":" + JsonOut::num(cxs[i]) + "," + JsonOut::str("y") + ":" + JsonOut::num(cys[i]) + "}";
        }
        J.body += "]";
        if (p.frame && p.frame_layers > 0 && p.frame_width > 0.0) {
            // Frame inner/outer bbox in bed coords (model-space ±(fid_outer + margin) and +frame_width).
            const double fid_outer        = fid_off + 2.5;
            const double frame_inner_edge = fid_outer + p.frame_margin;
            const double frame_outer_edge = frame_inner_edge + p.frame_width;
            const double frame_z_max      = 0.20 * static_cast<double>(p.frame_layers);
            meta << boost::format("; frame_enabled = 1\n");
            meta << boost::format("; frame_outer_bbox = X%1$.2f,Y%2$.2f X%3$.2f,Y%4$.2f\n")
                % (bed_center.x() - frame_outer_edge) % (bed_center.y() - frame_outer_edge)
                % (bed_center.x() + frame_outer_edge) % (bed_center.y() + frame_outer_edge);
            meta << boost::format("; frame_inner_bbox = X%1$.2f,Y%2$.2f X%3$.2f,Y%4$.2f\n")
                % (bed_center.x() - frame_inner_edge) % (bed_center.y() - frame_inner_edge)
                % (bed_center.x() + frame_inner_edge) % (bed_center.y() + frame_inner_edge);
            meta << boost::format("; frame_layers = %1%\n") % p.frame_layers;
            meta << boost::format("; frame_z_max_mm = %1$.2f\n") % frame_z_max;
            J.body += ",\n";
            J.w(JsonOut::str("frame") + ": {"
                + JsonOut::str("enabled") + ":true,"
                + JsonOut::str("outer_bbox") + ":{"
                    + JsonOut::str("x_min") + ":" + JsonOut::num(bed_center.x() - frame_outer_edge) + ","
                    + JsonOut::str("y_min") + ":" + JsonOut::num(bed_center.y() - frame_outer_edge) + ","
                    + JsonOut::str("x_max") + ":" + JsonOut::num(bed_center.x() + frame_outer_edge) + ","
                    + JsonOut::str("y_max") + ":" + JsonOut::num(bed_center.y() + frame_outer_edge) + "},"
                + JsonOut::str("inner_bbox") + ":{"
                    + JsonOut::str("x_min") + ":" + JsonOut::num(bed_center.x() - frame_inner_edge) + ","
                    + JsonOut::str("y_min") + ":" + JsonOut::num(bed_center.y() - frame_inner_edge) + ","
                    + JsonOut::str("x_max") + ":" + JsonOut::num(bed_center.x() + frame_inner_edge) + ","
                    + JsonOut::str("y_max") + ":" + JsonOut::num(bed_center.y() + frame_inner_edge) + "},"
                + JsonOut::str("layers") + ":" + std::to_string(p.frame_layers) + ","
                + JsonOut::str("z_max_mm") + ":" + JsonOut::num(frame_z_max)
                + "}");
        }
        J.body += "\n";
    }
    else if (tower_params && (type == CLICalibType::TempTower || type == CLICalibType::VolSpeedTower
                              || type == CLICalibType::PATower || type == CLICalibType::RetractionTower
                              || type == CLICalibType::VFATower)) {
        const auto &p = *tower_params;
        meta << boost::format("; tower_start = %1%\n") % p.start;
        meta << boost::format("; tower_end = %1%\n") % p.end;
        meta << boost::format("; tower_step = %1%\n") % p.step;
        auto blocks = compute_tower_blocks(type, p);
        const char* key = tower_value_key(type);
        for (size_t i = 0; i < blocks.size(); ++i) {
            const auto &b = blocks[i];
            meta << boost::format("; tower_block_%1% = z_low=%2$.2f, z_high=%3$.2f, %4%=%5%\n")
                % i % b.z_low % b.z_high % key % JsonOut::num(b.value);
        }

        // JSON
        J.w(JsonOut::str("params") + ": {"
            + JsonOut::str("start") + ":" + JsonOut::num(p.start) + ","
            + JsonOut::str("end")   + ":" + JsonOut::num(p.end)   + ","
            + JsonOut::str("step")  + ":" + JsonOut::num(p.step)  + "},\n");
        J.w(JsonOut::str("value_key") + ": " + JsonOut::str(key) + ",\n");
        J.w(JsonOut::str("blocks") + ": [\n"); J.indent = 2;
        for (size_t i = 0; i < blocks.size(); ++i) {
            J.w("{" + JsonOut::str("index") + ":" + std::to_string(i)
                + "," + JsonOut::str("z_low_mm") + ":" + JsonOut::num(blocks[i].z_low)
                + "," + JsonOut::str("z_high_mm") + ":" + JsonOut::num(blocks[i].z_high)
                + "," + JsonOut::str(key) + ":" + JsonOut::num(blocks[i].value) + "}");
            if (i + 1 < blocks.size()) J.body += ",";
            J.body += "\n";
        }
        J.indent = 1; J.w("]\n");
    }
    else if (flow_params && (type == CLICalibType::FlowRate_YOLO_Recommended
                             || type == CLICalibType::FlowRate_YOLO_Perfectionist
                             || type == CLICalibType::FlowRate_Pass1
                             || type == CLICalibType::FlowRate_Pass2)) {
        double cur_flowrate = 1.0;
        if (auto *fr = full_config.option<ConfigOptionFloatsNullable>("filament_flow_ratio"); fr && !fr->values.empty())
            cur_flowrate = fr->values[0];
        else if (auto *fr2 = full_config.option<ConfigOptionFloats>("filament_flow_ratio"); fr2 && !fr2->values.empty())
            cur_flowrate = fr2->values[0];
        auto blocks = compute_flow_blocks(model, cur_flowrate, flow_params->is_linear);
        meta << boost::format("; flow_seed = %1$.4f\n") % cur_flowrate;
        meta << boost::format("; is_linear = %1%\n") % (flow_params->is_linear ? 1 : 0);
        meta << boost::format("; block_count = %1%\n") % blocks.size();
        for (size_t i = 0; i < blocks.size(); ++i) {
            const auto &b = blocks[i];
            meta << boost::format("; flow_block_%1% = name=%2%, modifier=%3$.4f, print_flow_ratio=%4$.4f, X=%5$.2f, Y=%6$.2f\n")
                % i % b.name % b.modifier % b.print_flow_ratio % b.x % b.y;
        }

        // JSON
        J.w(JsonOut::str("flow_seed") + ": " + JsonOut::num(cur_flowrate) + ",\n");
        J.w(JsonOut::str("is_linear") + ": " + (flow_params->is_linear ? "true" : "false") + ",\n");
        J.w(JsonOut::str("blocks") + ": [\n"); J.indent = 2;
        for (size_t i = 0; i < blocks.size(); ++i) {
            J.w("{" + JsonOut::str("name") + ":" + JsonOut::str(blocks[i].name)
                + "," + JsonOut::str("modifier") + ":" + JsonOut::num(blocks[i].modifier)
                + "," + JsonOut::str("print_flow_ratio") + ":" + JsonOut::num(blocks[i].print_flow_ratio)
                + "," + JsonOut::str("x") + ":" + JsonOut::num(blocks[i].x)
                + "," + JsonOut::str("y") + ":" + JsonOut::num(blocks[i].y) + "}");
            if (i + 1 < blocks.size()) J.body += ",";
            J.body += "\n";
        }
        J.indent = 1; J.w("]\n");
    }
    meta << "; ---- END ORCA CALIBRATION METADATA ----\n";
    J.indent = 0; J.w("}\n");

    // ---- Inject metadata block into the gcode header ----
    {
        std::ifstream in(gcode_path);
        if (!in) return false;
        std::ostringstream body;
        body << in.rdbuf();
        in.close();
        std::string contents = body.str();
        const std::string marker = "; HEADER_BLOCK_START";
        size_t insert_pos = 0;
        size_t marker_pos = contents.find(marker);
        if (marker_pos != std::string::npos) {
            size_t line_end = contents.find('\n', marker_pos);
            insert_pos = (line_end != std::string::npos) ? line_end + 1 : marker_pos + marker.size();
        }
        contents.insert(insert_pos, meta.str());
        std::ofstream out(gcode_path, std::ios::trunc);
        if (!out) return false;
        out.write(contents.data(), contents.size());
    }

    // ---- Write JSON sidecar ----
    boost::filesystem::path json_path(gcode_path);
    json_path.replace_extension(".calib.json");
    {
        std::ofstream js(json_path.string(), std::ios::trunc);
        if (js) js << J.body;
    }
    BOOST_LOG_TRIVIAL(info) << "cli_emit_calib_outputs: wrote metadata into " << gcode_path
                            << " and sidecar " << json_path.string();
    return true;
}

}
