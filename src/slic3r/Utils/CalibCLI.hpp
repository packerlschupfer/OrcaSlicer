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
    ZOffsetPattern,
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

// ============================================================
// Z-offset / Live-Z calibration pattern (scanner-friendly).
// ============================================================
//
// Builds a procedural single-layer test plate optimized for AI-vision analysis of a flatbed
// scanner image. The output gcode contains:
//   - 4 fiducial marks at corner positions for auto-alignment / de-skew
//   - A 10mm scale bar with 1mm ticks for DPI cross-verification
//   - Zone S (solid): squish uniformity indicator (concentric infill pad)
//   - Zone G (gaps):  3 line-pairs at 0.5/0.6/0.8mm spacing — under-Z closes gaps, over-Z leaves all open
//   - Zone W (walls): 3 single-perimeter freestanding walls — sensitive to over-squish
//   - Zone C (corners): 4 small loops at the 4 plate corners — cross-bed Z uniformity
// Coordinates of every primitive are written to the gcode header as comments so downstream
// AI-vision tooling can map scan pixels to test zones without inferring from geometry.
//
// The plate is built by populating the (pre-loaded placeholder) Model with TriangleMesh box
// primitives — one ModelObject per primitive — with per-object configs that force each zone
// to print as intended (single perimeter for W, solid concentric for S, etc.). The CLI's
// standard --slice 0 pipeline handles PRINT_START, temperatures, and extrusion math; this
// generator only owns the geometry + per-object overrides + a single-layer hard cap via
// layer_height + initial_layer_print_height = 0.20mm.

struct CLIZCalParams {
    double plate_size = 100.0;     // mm, total plate footprint (square)
    double zone_size  = 30.0;      // mm, individual zone footprint for S/G/W
    bool   fiducials  = true;      // 4 corner fiducial marks for AI-vision auto-alignment
    bool   scale_bar  = true;      // 0-10mm scale bar with 1mm ticks for DPI verification
    bool   zone_labels = false;    // v1: skip stick-font labels (metadata comments suffice)
};

void cli_build_zcal_pattern(Model &model, DynamicPrintConfig &full_config, const CLIZCalParams &params);

// Post-process the sliced gcode to inject coordinate metadata comments. Called by OrcaSlicer.cpp
// after --slice 0 produces the gcode for the z-offset-pattern. Reads bed_center from the config
// to convert model-space coordinates to bed coordinates before writing. Returns true on success.
bool cli_inject_zcal_metadata(const std::string &gcode_path, const CLIZCalParams &params,
                              const DynamicPrintConfig &full_config);

}

#endif
