#ifndef slic3r_CalibCLI_hpp_
#define slic3r_CalibCLI_hpp_

#include <string>

#include "libslic3r/calib.hpp"
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
    TempTower,
    VolSpeedTower,
    PATower,
    RetractionTower,
    VFATower,
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
    double plate_size = 60.0;      // mm, total plate footprint (square) — v8 compact default
    double zone_size  = 18.0;      // mm, individual zone footprint for S/G/W — fits 3 zones in a row
    bool   fiducials  = true;      // 4 corner fiducial marks for AI-vision auto-alignment
    bool   scale_bar  = true;      // 0-10mm scale bar with 1mm ticks for DPI verification
    bool   zone_labels = false;    // v1: skip stick-font labels (metadata comments suffice)
    // ORCA: peelable structural frame outside the fiducial bbox. Calibration zones / fiducials /
    //       scale bar / corner-C loops stay strictly single-layer; only the frame is multi-layer.
    bool   frame = true;
    int    frame_layers = 2;       // 2 × 0.20mm = 0.40mm
    double frame_width = 3.0;
    double frame_margin = 5.0;
    // ORCA: connective struts + fragment IDs (struts brief, 2026-06-17). With the v8 grid mesh
    //       (see below), individual struts are largely redundant — the grid wires everything.
    //       Kept for backward-compat; default ON.
    bool   struts = true;
    bool   fiducial_ids = true;    // each fiducial becomes an L-shape with notch in a unique corner
    bool   zone_ids = true;        // 1/2/3 small dots beside S/G/W identifying which zone
    // ORCA v8: structural grid mesh covering the interior of the frame. Replaces single struts
    //         with a 4mm-pitch rectangular grid (horizontal + vertical lines, each 0.45mm wide
    //         single-perimeter). The grid bonds to every island via 0.5mm overlap at the
    //         boundaries — tension is distributed across many connections, no single failure point.
    //         Like rebar in concrete. Default ON.
    bool   grid = true;
    double grid_pitch_mm = 4.0;    // mm spacing between grid lines (horizontal + vertical)
    bool   grid_over_zones = false; // when false: lines stop at zone/fid/etc. borders (clean signal area)
                                    // when true:  lines cross zones (grid visible inside scan)
};

void cli_build_zcal_pattern(Model &model, DynamicPrintConfig &full_config, const CLIZCalParams &params);

// ============================================================
// Tower-shaped calibration tests (Temperature, Volumetric Speed, PA Tower).
// ============================================================
//
// All three use a pre-built tower mesh in resources/calib/ which the slicer modulates per-layer
// via Print::set_calib_params() — e.g. for the temp tower the gcode emitter inserts M104/M109
// commands at every 10mm Z transition; for the PA tower it inserts M572 (SET_PRESSURE_ADVANCE);
// for the vol-speed tower it ramps speed.
//
// The cli_apply_* functions handle the geometry side of what the GUI's Plater::calib_temp /
// calib_max_vol_speed / _calib_pa_tower do:
//   - load the resource mesh into the Model (the input-file pipeline does this)
//   - cut the tower at the top (and for temp, also the bottom) to match the user's range
//   - scale by nozzle ratio where the GUI does it
//   - apply per-object and print/filament config overrides line-for-line from the GUI
//
// Caller (OrcaSlicer.cpp) is responsible for calling Print::set_calib_params() with the result
// of cli_tower_get_calib_params() before slicing — that's what makes the test "active" at
// gcode-generation time.

struct CLITowerParams {
    double start = 0.0;  // tower start value (T_start °C / vs_start mm³/s / pa_start s)
    double end   = 0.0;  // tower end   value
    double step  = 0.0;  // step size between blocks (ignored for temp tower — fixed at 5°C)
};

void cli_apply_temp_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params);
void cli_apply_vol_speed_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params);
void cli_apply_pa_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params);
void cli_apply_retraction_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params);
void cli_apply_vfa_tower(Model &model, DynamicPrintConfig &full_config, const CLITowerParams &params);

// Translate the CLI tower params to the Calib_Params the Print engine consumes via
// set_calib_params(). For vol-speed tower this internally converts mm³/s → mm/s using the
// effective layer flow (mm³/mm). Returns true on success; out is left default on failure.
bool cli_tower_get_calib_params(CLICalibType type, const CLITowerParams &params,
                                const DynamicPrintConfig &full_config, Calib_Params &out);

// ============================================================
// Unified calibration outputs (gcode metadata + JSON sidecar).
// ============================================================
//
// After --slice 0 produces the gcode, this hook injects calibration metadata into the gcode
// header (machine-readable comments) AND writes a <basename>.calib.json file next to the gcode
// describing the test in fully structured form. Both are designed to be consumed by AI-vision
// analysis tooling that maps a scan back to test parameters without having to invert the
// calibration math from start/end/step.
//
// Coverage by calib type:
//   - Z-offset: per-zone coordinates + fiducial + scale-bar geometry (existing behavior).
//   - Tower types: per-block table mapping Z-range → block parameter
//     (temp °C / PA s / vol-speed mm³/s / retraction mm / VFA mm/s).
//   - Flow-rate (YOLO + Pass): per-block list of object name + modifier + resulting flow ratio.
//
// Returns true if the file was written; false if gcode_path is missing or unwritable.
bool cli_emit_calib_outputs(const std::string &gcode_path,
                            CLICalibType type,
                            const DynamicPrintConfig &full_config,
                            const CLIZCalParams *zcal_params,         // non-null only for ZOffsetPattern
                            const CLITowerParams *tower_params,       // non-null only for tower types
                            const CLIFlowRateParams *flow_params,     // non-null only for flow types
                            const Model *model);                       // for flow per-block listing

}

#endif
