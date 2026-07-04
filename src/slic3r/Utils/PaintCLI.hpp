// PaintCLI.hpp — CLI paint-inspection primitives.
//
// Backs the --inspect-paint CLI action. Reads the per-facet enforcer /
// blocker / extruder / fuzzy-skin state that OrcaSlicer stores on every
// ModelVolume (supports, seam, MMU color, fuzzy-skin) and emits a
// structured JSON summary — facet count, surface area, and mesh-local
// bounding box per state — so CI / scripted / AI tooling can reason
// about existing paint on a .3mf without opening the GUI.
//
// Coordinates are mesh-local (each volume's own frame), matching the
// frame that the paint gizmos operate in.
#ifndef slic3r_PaintCLI_hpp_
#define slic3r_PaintCLI_hpp_

#include <iosfwd>
#include <string>
#include <vector>

namespace Slic3r {
class Model;

namespace PaintCLI {

// `source_paths` lists every input file; the CLI merges them into one Model.
void inspect_to_json(const Model &model, const std::vector<std::string> &source_paths,
                     std::ostream &out);

// Apply a --paint-supports spec.json to `model` in place. The spec selects
// one object + optional volume, then applies an ordered list of regions;
// each region is (kind, predicate) where kind is enforcer/blocker/clear and
// the predicate filters source triangles by bbox / z_min / z_max /
// normal_below (n·(-Z) ≥ threshold, i.e. facing down). Emits a JSON
// coverage report to `report_out`. Returns true if every region matched
// ≥1 facet; false if any matched 0 (the caller can use this together with
// --strict to fail the whole action).
//
// All coordinates in the spec are mesh-local. Frame conversion is out of
// scope for this step; a future --frame=world flag will add it.
bool apply_supports_spec(Model &model, const std::string &spec_path,
                         const std::string &source_path,
                         std::ostream &report_out);

} // namespace PaintCLI
} // namespace Slic3r

#endif
