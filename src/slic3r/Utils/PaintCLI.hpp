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

} // namespace PaintCLI
} // namespace Slic3r

#endif
