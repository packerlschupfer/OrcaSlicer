// PaintCLI.hpp — CLI paint-inspection primitives.
//
// Reads the per-facet enforcer/blocker bitstreams that OrcaSlicer stores on
// every ModelVolume (supports, seam, MMU color, fuzzy-skin) and emits a
// structured JSON summary — facet count, surface area, and mesh-local
// bounding box per state — for AI/CLI tooling to reason about existing
// paint without loading the GUI.
//
// Coordinates are mesh-local (each volume's own frame), matching the frame
// the paint gizmos operate in and the frame --ground-face-point uses.
#ifndef slic3r_PaintCLI_hpp_
#define slic3r_PaintCLI_hpp_

#include <iosfwd>
#include <string>

namespace Slic3r {
class Model;

namespace PaintCLI {

void inspect_to_json(const Model &model, const std::string &source_path,
                     std::ostream &out);

// Apply a --paint spec.json to `model` in place. The spec selects one
// object + optional volume and a paint layer (supports/seam/mmu/fuzzy,
// default supports), then applies an ordered list of regions. Each
// region is (kind, predicate); predicates filter source triangles by
// bbox / z_min / z_max / normal_below, plus optional reachable_from
// (point-seeded flood-fill) and connected_component (keep only the
// largest connected match subset). Coordinates are mesh-local by
// default; set "frame": "world" for bed coords. Emits a JSON coverage
// report to `report_out`. Returns true if every region matched ≥1
// facet; false if any matched 0 (the caller can pair this with
// --strict to fail the action).
bool apply_spec(Model &model, const std::string &spec_path,
                const std::string &source_path,
                std::ostream &report_out);

// Options for render_paint. All optional; sensible defaults for a
// support-focused headless verify.
struct RenderPaintOpts {
    std::string layer          = "supports";  // supports | seam | mmu | fuzzy
    std::string view           = "top";       // top (–Z) | bottom (+Z)
    int         width_px       = 800;
    int         height_px      = 800;
    // Optional horizontal-cut clip. If clip_z_below is finite, only paint
    // at world_z ≤ clip_z_below is rendered; the ray advances past hidden
    // material. clip_z_above is the mirror (paint at world_z ≥ Z).
    double      clip_z_below   = std::numeric_limits<double>::infinity();
    double      clip_z_above   = -std::numeric_limits<double>::infinity();
};

// Emit a PPM (P6 binary) top-down (or bottom-up) view of the loaded model
// with painted regions color-coded. Sub-facet paint (post-split state) is
// resolved by walking the TriangleSelector split tree. Ray misses are
// black; unpainted mesh is gray-shaded by Lambert dot(normal, -view_dir);
// enforcer is red, blocker is blue, MMU extruders map to a palette,
// fuzzy_skin is green. Header line "P6\nW H\n255\n" then RGB bytes.
//
// Returns true on success.
bool render_paint(const Model &model, const std::string &out_path,
                  const RenderPaintOpts &opts, std::ostream &report_out);

} // namespace PaintCLI
} // namespace Slic3r

#endif
