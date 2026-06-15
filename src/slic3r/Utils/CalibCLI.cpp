#include "CalibCLI.hpp"

#include <algorithm>
#include <cmath>

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

        //ORCA: brim toggle from PR #13548 — when enabled, switch each block-pair object to an outer
        //      brim with the requested width and zero object-gap (matching the dialog's behavior).
        if (params.brim_enabled) {
            mo->config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btOuterOnly));
            mo->config.set_key_value("brim_width", new ConfigOptionFloat(params.brim_width));
            mo->config.set_key_value("brim_object_gap", new ConfigOptionFloat(0.0));
        }
    }

    //ORCA: printer-level override mirrors the GUI implementation.
    full_config.set_key_value("resonance_avoidance", new ConfigOptionBool(false));

    BOOST_LOG_TRIVIAL(info) << boost::format("cli_apply_flowrate_calib: pass=%1% is_linear=%2% pattern=%3% brim=%4% objects=%5%")
        % params.pass % params.is_linear % static_cast<int>(params.pattern) % params.brim_enabled % model.objects.size();
}

}
