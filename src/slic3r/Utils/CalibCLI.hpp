#ifndef slic3r_CalibCLI_hpp_
#define slic3r_CalibCLI_hpp_

#include <string>

#include "libslic3r/Model.hpp"
#include "libslic3r/PrintConfig.hpp"


namespace Slic3r {

// ORCA: CLI calibration support. These free functions are GUI-independent equivalents of the
//       Plater::calib_* methods used by the Calibration Wizard. They take a Model (already
//       populated from the calibration .3mf resource), a merged DynamicPrintConfig (the CLI's
//       m_print_config after --load-settings / --load-filaments processing), and calibration-
//       specific parameters. After calling, the Model's per-object configs and the passed
//       DynamicPrintConfig contain the same overrides the GUI wizard would have applied. The
//       CLI then continues into the standard --slice pipeline, producing G-code byte-identical
//       (modulo timestamp comments) to a GUI wizard run with the same selections.
//
//       v1 scope: flow-rate tests. Temperature tower / pressure advance tower follow once the
//       flow path is validated against acceptance criterion #1 (CLI vs GUI byte-identical diff).

enum class CLICalibType {
    NoCalib,
    FlowRate_YOLO_Recommended,
    FlowRate_YOLO_Perfectionist,
    FlowRate_Pass1,
    FlowRate_Pass2,
};

// Maps a CLI string ("flow-yolo-recommended", etc.) to the enum + a 3MF resource path.
// Returns CLICalibType::None on unknown string.
CLICalibType cli_calib_type_from_string(const std::string &name);

// Resolve the model resource path inside resources_dir() for a given calibration type.
// Returns empty string for None.
std::string cli_calib_resource_path(CLICalibType type);

// Flow rate (YOLO and Pass1/Pass2) — apply per-object config overrides + geometry scale.
struct CLIFlowRateParams {
    int            pass = 1;                     // 1 = coarse / YOLO-recommended, 2 = fine / YOLO-perfectionist
    bool           is_linear = true;              // true = YOLO variants, false = Pass1/Pass2 variants
    InfillPattern  pattern = ipArchimedeanChords; // top-surface pattern
    bool           brim_enabled = false;
    double         brim_width = 2.0;              // mm, only used when brim_enabled
    double         brim_extra_gap = 0.0;          // mm, only used when brim_enabled
    int            max_blocks = -1;               // -1 = all; otherwise keep N blocks closest to modifier=0
    double         max_modifier = -1.0;           // <=0 = no filter; otherwise keep blocks with |modifier| <= this
};
void cli_apply_flowrate_calib(Model &model, DynamicPrintConfig &full_config, const CLIFlowRateParams &params);

// Returns (pass, is_linear) for a CLI calib type. Caller fills in pattern/brim from CLI flags.
void cli_flowrate_params_for_type(CLICalibType type, int &pass, bool &is_linear);

// Map the CLI's --flow-pattern string to an InfillPattern. Default ipArchimedeanChords on bad input.
InfillPattern cli_parse_flow_pattern(const std::string &name);

}

#endif
