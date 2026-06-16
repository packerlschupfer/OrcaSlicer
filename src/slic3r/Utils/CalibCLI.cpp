#include "CalibCLI.hpp"

#include <algorithm>
#include <cmath>
#include <clocale>
#include <cstdlib>
#include <string>

#include <boost/log/trivial.hpp>
#include <boost/format.hpp>

#include "libslic3r/Flow.hpp"
#include "libslic3r/Utils.hpp"

namespace Slic3r {

CLICalibType cli_calib_type_from_string(const std::string &name)
{
    if (name == "flow-yolo-recommended")    return CLICalibType::FlowRate_YOLO_Recommended;
    if (name == "flow-yolo-perfectionist")  return CLICalibType::FlowRate_YOLO_Perfectionist;
    if (name == "flow-pass1")               return CLICalibType::FlowRate_Pass1;
    if (name == "flow-pass2")               return CLICalibType::FlowRate_Pass2;
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

}
