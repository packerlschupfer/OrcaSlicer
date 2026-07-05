// PaintCLI.cpp — CLI paint-inspection primitives. See PaintCLI.hpp.
#include "PaintCLI.hpp"

#include "libslic3r/AABBTreeIndirect.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/TriangleSelector.hpp"

#include <boost/log/trivial.hpp>
#include <boost/nowide/fstream.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r {
namespace PaintCLI {

namespace {

using json = nlohmann::json;

double its_surface_area(const indexed_triangle_set &its)
{
    double total = 0.0;
    for (const stl_triangle_vertex_indices &t : its.indices) {
        const Vec3f &a = its.vertices[t(0)];
        const Vec3f &b = its.vertices[t(1)];
        const Vec3f &c = its.vertices[t(2)];
        total += 0.5 * (b - a).cross(c - a).norm();
    }
    return total;
}

// Bbox over triangle-referenced vertices only. get_facets_strict() returns
// an itset with the full source vertex list — using bounding_box() on it
// would report the whole mesh's bbox even when only a few facets are painted.
BoundingBoxf3 its_referenced_bbox(const indexed_triangle_set &its)
{
    BoundingBoxf3 bb;
    bool first = true;
    for (const stl_triangle_vertex_indices &t : its.indices) {
        for (int k = 0; k < 3; ++k) {
            const Vec3d v = its.vertices[t(k)].cast<double>();
            if (first) { bb.min = bb.max = v; first = false; }
            else       bb.merge(v);
        }
    }
    return bb;
}

json vec3_to_json(const Vec3d &v)
{
    return json::array({ v.x(), v.y(), v.z() });
}

json bbox_to_json(const BoundingBoxf3 &bb)
{
    return {
        { "min",  vec3_to_json(bb.min) },
        { "max",  vec3_to_json(bb.max) },
        { "size", vec3_to_json(Vec3d(bb.max - bb.min)) },
    };
}

// One (layer, state) row — empty ones are omitted at the caller level.
json state_entry(const std::string &label, const indexed_triangle_set &its)
{
    return {
        { "state",    label },
        { "facets",   its.indices.size() },
        { "area_mm2", its_surface_area(its) },
        { "bbox",     bbox_to_json(its_referenced_bbox(its)) },
    };
}

// Iterate the states relevant to one FacetsAnnotation kind, collecting
// non-empty entries. Empty layer → {"empty": true}. `n_facets_out` is the
// running total of painted facets — bumped for the summary.
json inspect_layer(const ModelVolume &mv, const FacetsAnnotation &fa,
                   const std::vector<std::pair<EnforcerBlockerType, std::string>> &states,
                   size_t &n_facets_out)
{
    if (fa.empty())
        return { { "empty", true } };

    json entries = json::array();
    for (const auto &st : states) {
        if (!fa.has_facets(mv, st.first))
            continue;
        indexed_triangle_set its = fa.get_facets_strict(mv, st.first);
        if (its.indices.empty())
            continue;
        n_facets_out += its.indices.size();
        entries.push_back(state_entry(st.second, its));
    }
    return {
        { "empty",  entries.empty() },
        { "states", std::move(entries) },
    };
}

const std::vector<std::pair<EnforcerBlockerType, std::string>> &supports_states()
{
    static const std::vector<std::pair<EnforcerBlockerType, std::string>> s = {
        { EnforcerBlockerType::ENFORCER, "ENFORCER" },
        { EnforcerBlockerType::BLOCKER,  "BLOCKER"  },
    };
    return s;
}

const std::vector<std::pair<EnforcerBlockerType, std::string>> &fuzzy_states()
{
    // FUZZY_SKIN is an enum alias for ENFORCER; the layer is single-state.
    static const std::vector<std::pair<EnforcerBlockerType, std::string>> s = {
        { EnforcerBlockerType::FUZZY_SKIN, "FUZZY_SKIN" },
    };
    return s;
}

const std::vector<std::pair<EnforcerBlockerType, std::string>> &mmu_states()
{
    static std::vector<std::pair<EnforcerBlockerType, std::string>> s = []{
        std::vector<std::pair<EnforcerBlockerType, std::string>> v;
        for (int i = 1; i <= int(EnforcerBlockerType::ExtruderMax); ++i)
            v.emplace_back(EnforcerBlockerType(i), "extruder_" + std::to_string(i));
        return v;
    }();
    return s;
}

} // namespace

void inspect_to_json(const Model &model, const std::string &source_path,
                     std::ostream &out)
{
    json root;
    root["source"] = source_path;
    root["frame"]  = "mesh_local";
    root["note"]   = "Coordinates are mesh-local (each volume's own frame). "
                     "Paint gizmos and the future --paint predicates "
                     "operate in this frame.";

    json objects = json::array();
    size_t total_objects = 0, total_volumes = 0, total_painted = 0, total_facets = 0;

    for (size_t oi = 0; oi < model.objects.size(); ++oi) {
        const ModelObject *mo = model.objects[oi];
        if (!mo) continue;
        ++total_objects;

        json obj;
        obj["index"] = oi;
        obj["name"]  = mo->name;

        json volumes = json::array();
        for (size_t vi = 0; vi < mo->volumes.size(); ++vi) {
            const ModelVolume *mv = mo->volumes[vi];
            if (!mv) continue;
            ++total_volumes;

            const indexed_triangle_set &its = mv->mesh().its;
            json vol;
            vol["index"]    = vi;
            vol["name"]     = mv->name;
            vol["n_facets"] = its.indices.size();
            vol["is_model_part"] = mv->is_model_part();
            vol["bbox_mesh_local"] = bbox_to_json(bounding_box(its));

            size_t vol_painted = 0;
            json paints;
            paints["supports"]         = inspect_layer(*mv, mv->supported_facets,
                                                      supports_states(),  vol_painted);
            paints["seam"]             = inspect_layer(*mv, mv->seam_facets,
                                                      supports_states(),  vol_painted);
            paints["mmu_segmentation"] = inspect_layer(*mv, mv->mmu_segmentation_facets,
                                                      mmu_states(),       vol_painted);
            paints["fuzzy_skin"]       = inspect_layer(*mv, mv->fuzzy_skin_facets,
                                                      fuzzy_states(),     vol_painted);
            vol["paints"] = std::move(paints);
            vol["painted_facets_total"] = vol_painted;

            if (vol_painted > 0) ++total_painted;
            total_facets += vol_painted;

            volumes.push_back(std::move(vol));
        }
        obj["volumes"] = std::move(volumes);
        objects.push_back(std::move(obj));
    }
    root["objects"] = std::move(objects);
    root["summary"] = {
        { "objects",                    total_objects },
        { "volumes",                    total_volumes },
        { "volumes_with_paint",         total_painted },
        { "painted_facets_total",       total_facets  },
    };

    out << root.dump(2) << std::endl;
}

// =============================== --paint ===========================

namespace {

// One predicate per region. All present fields AND together; a field
// left unset is a no-op.
struct Predicate {
    std::optional<BoundingBoxf3> bbox;
    std::optional<double>        z_min;
    std::optional<double>        z_max;
    std::optional<double>        normal_below;   // dot(n, -Z) ≥ this  ⇒ facing down
    std::optional<Vec3d>         reachable_from; // flood-fill from nearest source triangle;
                                                 // only neighbors that also match the other
                                                 // fields propagate. Isolates one connected
                                                 // face when combined with normal/bbox.
    bool                         connected_component = false; // after matching, keep only
                                                              // the largest connected subset
    bool                         all = false;    // shortcut: match every facet
};

// 2D shape projected onto the mesh — mirrors the GUI's line/rect/polygon
// painter tools. When present, the shape replaces the whole-triangle
// predicate: samples inside the shape are raycast onto the mesh and
// splatted with a Circle cursor, which subdivides triangles at the
// cursor edge (so the painted edge follows the shape, not the facet).
enum class ProjectAxis { NegZ, PosZ };  // from above (default) or from below

struct Shape {
    enum class Kind { None, Line, Rect, Polygon };
    Kind kind = Kind::None;

    // For all kinds: 2D points in the projection plane (XY when project=NegZ/PosZ).
    // Line uses points[0..1] + width_mm.
    // Rect uses points[0..1] as min/max corners.
    // Polygon uses points[0..n-1] as vertices (closed implicitly).
    std::vector<Vec2d> points;
    double             width_mm         = 0.0;   // line half-width scaled to 2× for rendering
    ProjectAxis        project          = ProjectAxis::NegZ;
    double             sample_spacing_mm = 0.5;  // grid density
    double             cursor_radius_mm = 0.0;   // splat radius (default: sample_spacing * 0.9)
};

struct Region {
    std::string         kind_str; // for report
    EnforcerBlockerType state;    // ENFORCER / BLOCKER / NONE (=clear)
    Predicate           pred;
    Shape               shape;    // if kind != None, splat path is used instead
                                  // of the whole-triangle predicate
};

Predicate parse_predicate(const json &sel)
{
    Predicate p;
    if (sel.value("all", false)) { p.all = true; return p; }
    if (sel.contains("bbox")) {
        const auto &b = sel["bbox"];
        BoundingBoxf3 bb;
        bb.min = Vec3d(b["min"][0].get<double>(), b["min"][1].get<double>(), b["min"][2].get<double>());
        bb.max = Vec3d(b["max"][0].get<double>(), b["max"][1].get<double>(), b["max"][2].get<double>());
        p.bbox = bb;
    }
    if (sel.contains("z_min"))        p.z_min        = sel["z_min"].get<double>();
    if (sel.contains("z_max"))        p.z_max        = sel["z_max"].get<double>();
    if (sel.contains("normal_below")) p.normal_below = sel["normal_below"].get<double>();
    if (sel.contains("reachable_from")) {
        const auto &r = sel["reachable_from"];
        p.reachable_from = Vec3d(r[0].get<double>(), r[1].get<double>(), r[2].get<double>());
    }
    if (sel.value("connected_component", false))
        p.connected_component = true;
    return p;
}

ProjectAxis parse_project(const std::string &s)
{
    if (s.empty() || s == "-Z") return ProjectAxis::NegZ;
    if (s == "+Z")               return ProjectAxis::PosZ;
    throw std::runtime_error("--paint: unsupported project axis \"" + s +
                             "\" (only -Z / +Z supported in this build)");
}

Shape parse_shape(const json &sel)
{
    Shape s;
    if (sel.contains("project"))          s.project = parse_project(sel["project"].get<std::string>());
    if (sel.contains("sample_spacing_mm")) s.sample_spacing_mm = sel["sample_spacing_mm"].get<double>();
    if (sel.contains("cursor_radius_mm"))  s.cursor_radius_mm  = sel["cursor_radius_mm"].get<double>();

    if (sel.contains("line")) {
        const auto &L = sel["line"];
        s.kind = Shape::Kind::Line;
        s.points = { Vec2d(L["a"][0].get<double>(), L["a"][1].get<double>()),
                     Vec2d(L["b"][0].get<double>(), L["b"][1].get<double>()) };
        s.width_mm = L.value("width_mm", 1.0);
    } else if (sel.contains("rect")) {
        const auto &R = sel["rect"];
        s.kind = Shape::Kind::Rect;
        s.points = { Vec2d(R["min"][0].get<double>(), R["min"][1].get<double>()),
                     Vec2d(R["max"][0].get<double>(), R["max"][1].get<double>()) };
    } else if (sel.contains("polygon")) {
        const auto &P = sel["polygon"];
        if (!P.is_array() || P.size() < 3)
            throw std::runtime_error("--paint: polygon needs ≥ 3 vertices");
        s.kind = Shape::Kind::Polygon;
        s.points.reserve(P.size());
        for (const auto &p : P) s.points.emplace_back(p[0].get<double>(), p[1].get<double>());
    }
    if (s.cursor_radius_mm == 0.0) s.cursor_radius_mm = s.sample_spacing_mm * 0.9;
    return s;
}

bool predicate_matches(const Predicate &p, const Vec3d &centroid, const Vec3d &normal)
{
    if (p.all)                                                        return true;
    if (p.bbox        && !p.bbox->contains(centroid))                 return false;
    if (p.z_min       && centroid.z() < *p.z_min)                     return false;
    if (p.z_max       && centroid.z() > *p.z_max)                     return false;
    if (p.normal_below && (-normal.z()) < *p.normal_below)            return false;
    return true;
}

// Per-triangle 3-neighbor array, indexed by shared edge. Value -1 for a
// boundary edge with no adjacent triangle. Built once per volume and
// reused across all regions of a spec.
using FaceNeighbors = std::vector<std::array<int, 3>>;

FaceNeighbors compute_face_neighbors(const indexed_triangle_set &its)
{
    // Edge key = (min_vertex, max_vertex) — canonical ordering so
    // both incident triangles hash into the same slot.
    using EdgeKey = std::pair<int, int>;
    std::map<EdgeKey, std::array<int, 2>> edges; // -> up to two owner triangles

    auto canonical = [](int u, int v) -> EdgeKey {
        return u < v ? EdgeKey{u, v} : EdgeKey{v, u};
    };

    for (size_t i = 0; i < its.indices.size(); ++i) {
        const auto &t = its.indices[i];
        for (int k = 0; k < 3; ++k) {
            EdgeKey key = canonical(t(k), t((k + 1) % 3));
            auto it = edges.find(key);
            if (it == edges.end()) edges[key] = { int(i), -1 };
            else if (it->second[1] == -1) it->second[1] = int(i);
            // Non-manifold edges (>2 triangles) — silently keep the first two.
        }
    }

    FaceNeighbors neighbors(its.indices.size(), {-1, -1, -1});
    for (size_t i = 0; i < its.indices.size(); ++i) {
        const auto &t = its.indices[i];
        for (int k = 0; k < 3; ++k) {
            EdgeKey key = canonical(t(k), t((k + 1) % 3));
            const auto &owners = edges[key];
            int other = owners[0] == int(i) ? owners[1] : owners[0];
            neighbors[i][k] = other;
        }
    }
    return neighbors;
}

// Find source triangle whose centroid is closest to `point`. Used to seed
// flood fill; caller decides whether that seed also satisfies the region's
// other predicates.
int find_seed_triangle(const indexed_triangle_set &its, const Vec3d &point)
{
    int    best     = -1;
    double best_d2  = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < its.indices.size(); ++i) {
        const auto &t = its.indices[i];
        const Vec3d centroid = ((its.vertices[t(0)] + its.vertices[t(1)] + its.vertices[t(2)]).cast<double>()) / 3.0;
        const double d2 = (centroid - point).squaredNorm();
        if (d2 < best_d2) { best_d2 = d2; best = int(i); }
    }
    return best;
}

// Return the size and members of the largest connected component of
// `matched`. Two matched triangles are connected if they share an edge.
std::vector<bool> largest_connected_component(const std::vector<bool> &matched,
                                              const FaceNeighbors     &neighbors)
{
    const int n = int(matched.size());
    std::vector<bool> best(n, false);
    std::vector<bool> visited(n, false);
    int best_size = 0;
    std::vector<int> component;
    for (int seed = 0; seed < n; ++seed) {
        if (!matched[seed] || visited[seed]) continue;
        component.clear();
        std::queue<int> q; q.push(seed); visited[seed] = true;
        while (!q.empty()) {
            int t = q.front(); q.pop();
            component.push_back(t);
            for (int nb : neighbors[t])
                if (nb >= 0 && matched[nb] && !visited[nb]) { visited[nb] = true; q.push(nb); }
        }
        if (int(component.size()) > best_size) {
            best_size = int(component.size());
            std::fill(best.begin(), best.end(), false);
            for (int t : component) best[t] = true;
        }
    }
    return best;
}

enum class Layer { Supports, Seam, Mmu, Fuzzy };

Layer parse_layer(const std::string &s)
{
    if (s.empty() || s == "supports") return Layer::Supports;
    if (s == "seam")                  return Layer::Seam;
    if (s == "mmu")                   return Layer::Mmu;
    if (s == "fuzzy" || s == "fuzzy_skin") return Layer::Fuzzy;
    throw std::runtime_error("--paint: unknown layer \"" + s +
                             "\" (must be supports|seam|mmu|fuzzy)");
}

const char *layer_str(Layer l)
{
    switch (l) {
    case Layer::Supports: return "supports";
    case Layer::Seam:     return "seam";
    case Layer::Mmu:      return "mmu";
    case Layer::Fuzzy:    return "fuzzy";
    }
    return "?";
}

FacetsAnnotation &layer_facets(ModelVolume &mv, Layer l)
{
    switch (l) {
    case Layer::Supports: return mv.supported_facets;
    case Layer::Seam:     return mv.seam_facets;
    case Layer::Mmu:      return mv.mmu_segmentation_facets;
    case Layer::Fuzzy:    return mv.fuzzy_skin_facets;
    }
    return mv.supported_facets; // unreachable
}

// ────────────────────── 2D shape sampling helpers ──────────────────────

BoundingBoxf shape_bbox_2d(const Shape &s)
{
    BoundingBoxf b;
    if (s.points.empty()) return b;
    b.min = b.max = s.points.front();
    for (const Vec2d &p : s.points) {
        b.min = b.min.cwiseMin(p);
        b.max = b.max.cwiseMax(p);
    }
    if (s.kind == Shape::Kind::Line) {
        // pad by width_mm on both sides
        b.min -= Vec2d(s.width_mm, s.width_mm);
        b.max += Vec2d(s.width_mm, s.width_mm);
    }
    return b;
}

bool point_in_polygon_2d(const std::vector<Vec2d> &poly, const Vec2d &p)
{
    // Even-odd rule.
    bool inside = false;
    const size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const Vec2d &a = poly[i], &b = poly[j];
        if (((a.y() > p.y()) != (b.y() > p.y())) &&
            (p.x() < (b.x() - a.x()) * (p.y() - a.y()) / (b.y() - a.y() + 1e-30) + a.x()))
            inside = !inside;
    }
    return inside;
}

double dist2_point_to_segment_2d(const Vec2d &p, const Vec2d &a, const Vec2d &b)
{
    Vec2d ab = b - a;
    double len2 = ab.squaredNorm();
    if (len2 < 1e-20) return (p - a).squaredNorm();
    double t = std::clamp((p - a).dot(ab) / len2, 0.0, 1.0);
    Vec2d proj = a + t * ab;
    return (p - proj).squaredNorm();
}

bool point_in_shape_2d(const Shape &s, const Vec2d &p)
{
    switch (s.kind) {
    case Shape::Kind::Rect:
        return p.x() >= s.points[0].x() && p.x() <= s.points[1].x() &&
               p.y() >= s.points[0].y() && p.y() <= s.points[1].y();
    case Shape::Kind::Line: {
        const double r = s.width_mm * 0.5;
        return dist2_point_to_segment_2d(p, s.points[0], s.points[1]) <= r * r;
    }
    case Shape::Kind::Polygon:
        return point_in_polygon_2d(s.points, p);
    case Shape::Kind::None:
        return false;
    }
    return false;
}

EnforcerBlockerType parse_kind(const std::string &s, Layer layer)
{
    if (s == "clear") return EnforcerBlockerType::NONE;
    switch (layer) {
    case Layer::Supports:
    case Layer::Seam:
        if (s == "enforcer") return EnforcerBlockerType::ENFORCER;
        if (s == "blocker")  return EnforcerBlockerType::BLOCKER;
        throw std::runtime_error("--paint: kind \"" + s + "\" is not valid for layer \"" +
                                 layer_str(layer) + "\" (must be enforcer|blocker|clear)");
    case Layer::Fuzzy:
        if (s == "fuzzy_skin" || s == "enforcer") return EnforcerBlockerType::FUZZY_SKIN;
        throw std::runtime_error("--paint: kind \"" + s + "\" is not valid for layer "
                                 "\"fuzzy\" (must be fuzzy_skin|clear)");
    case Layer::Mmu: {
        if (s.rfind("extruder_", 0) == 0) {
            const int n = std::atoi(s.c_str() + 9);
            if (n >= 1 && n <= int(EnforcerBlockerType::ExtruderMax))
                return EnforcerBlockerType(n);
        }
        throw std::runtime_error("--paint: kind \"" + s + "\" is not valid for layer "
                                 "\"mmu\" (must be extruder_1..extruder_16 or clear)");
    }
    }
    return EnforcerBlockerType::NONE;
}

// Locate the target volume(s). Returns pointers into the model; does not
// take ownership. Empty result → the caller reports a warning.
std::vector<ModelVolume *> select_volumes(Model &model, const json &spec)
{
    std::vector<ModelVolume *> out;

    // Object selector: {index: N} or {name: "…"}. Default = first object.
    ModelObject *chosen_obj = nullptr;
    const json &osel = spec.value("object", json::object());
    if (osel.contains("index")) {
        int idx = osel["index"].get<int>();
        if (idx >= 0 && idx < int(model.objects.size())) chosen_obj = model.objects[idx];
    } else if (osel.contains("name")) {
        const std::string want = osel["name"].get<std::string>();
        for (ModelObject *o : model.objects) if (o && o->name == want) { chosen_obj = o; break; }
    } else if (!model.objects.empty()) {
        chosen_obj = model.objects.front();
    }
    if (!chosen_obj) return out;

    // Volume selector: {index: N} or {name: "…"} or none = every is_model_part.
    if (spec.contains("volume")) {
        const json &vsel = spec["volume"];
        if (vsel.contains("index")) {
            int idx = vsel["index"].get<int>();
            if (idx >= 0 && idx < int(chosen_obj->volumes.size()))
                out.push_back(chosen_obj->volumes[idx]);
        } else if (vsel.contains("name")) {
            const std::string want = vsel["name"].get<std::string>();
            for (ModelVolume *v : chosen_obj->volumes)
                if (v && v->name == want) { out.push_back(v); break; }
        }
    } else {
        for (ModelVolume *v : chosen_obj->volumes)
            if (v && v->is_model_part()) out.push_back(v);
    }
    return out;
}

} // namespace

bool apply_spec(Model &model, const std::string &spec_path,
                         const std::string &source_path,
                         std::ostream &report_out)
{
    json spec;
    {
        boost::nowide::ifstream in(spec_path.c_str());
        if (!in) {
            BOOST_LOG_TRIVIAL(error) << "--paint: cannot open " << spec_path;
            report_out << json({{ "error", "cannot_open_spec" }, { "path", spec_path }}).dump(2) << std::endl;
            return false;
        }
        try { in >> spec; }
        catch (const std::exception &e) {
            BOOST_LOG_TRIVIAL(error) << "--paint: invalid JSON in " << spec_path << ": " << e.what();
            report_out << json({{ "error", "invalid_json" }, { "path", spec_path }, { "detail", e.what() }}).dump(2) << std::endl;
            return false;
        }
    }

    // Layer + frame — both optional, defaulting to supports / mesh_local for
    // backward compatibility with the pre-step-4 spec shape.
    Layer layer;
    try { layer = parse_layer(spec.value("layer", std::string{"supports"})); }
    catch (const std::exception &e) {
        report_out << json({{ "error", "invalid_layer" }, { "detail", e.what() }}).dump(2) << std::endl;
        return false;
    }
    const std::string frame = spec.value("frame", std::string{"mesh_local"});
    if (frame != "mesh_local" && frame != "world") {
        report_out << json({{ "error", "invalid_frame" }, { "detail", "frame must be mesh_local|world" }}).dump(2) << std::endl;
        return false;
    }
    const bool frame_is_world = (frame == "world");

    // Optional clip plane — mirrors the GUI section-view. A facet whose
    // sampled hit point lies on the "clipped" side is discarded from the
    // shape splat / inspect pass. Two shorthands accepted:
    //   { "z_below": Z }  →  clip everything with world_z > Z  (paint only ≤ Z)
    //   { "z_above": Z }  →  clip everything with world_z < Z  (paint only ≥ Z)
    // Or the general form { "plane_normal": [x,y,z], "offset": d } where the
    // clipped side is n·p > d.
    TriangleSelector::ClippingPlane clip{}; // inactive by default
    if (spec.contains("clip")) {
        const json &c = spec["clip"];
        if (c.contains("z_below")) {
            clip.normal = { 0.f, 0.f, 1.f };
            clip.offset = float(c["z_below"].get<double>());
        } else if (c.contains("z_above")) {
            clip.normal = { 0.f, 0.f, -1.f };
            clip.offset = float(-c["z_above"].get<double>());
        } else if (c.contains("plane_normal")) {
            const auto &n = c["plane_normal"];
            clip.normal = { float(n[0].get<double>()), float(n[1].get<double>()), float(n[2].get<double>()) };
            clip.offset = float(c.value("offset", 0.0));
        }
    }

    // Parse regions upfront so errors surface before we touch the model.
    std::vector<Region> regions;
    if (!spec.contains("regions") || !spec["regions"].is_array()) {
        report_out << json({{ "error", "regions_missing_or_not_array" }}).dump(2) << std::endl;
        return false;
    }
    try {
        for (const json &r : spec["regions"]) {
            Region reg;
            reg.kind_str = r.value("kind", "");
            reg.state    = parse_kind(reg.kind_str, layer);
            reg.pred     = parse_predicate(r.value("select", json::object()));
            reg.shape    = parse_shape(r.value("select", json::object()));
            regions.push_back(std::move(reg));
        }
    } catch (const std::exception &e) {
        report_out << json({{ "error", "invalid_region" }, { "detail", e.what() }}).dump(2) << std::endl;
        return false;
    }

    std::vector<ModelVolume *> targets = select_volumes(model, spec);
    if (targets.empty()) {
        report_out << json({
            { "source",   source_path },
            { "spec",     spec_path   },
            { "error",    "no_volume_matched" },
            { "hint",     "check object/volume selector; --inspect-paint lists valid names/indices" },
        }).dump(2) << std::endl;
        return false;
    }

    json report;
    report["source"] = source_path;
    report["spec"]   = spec_path;
    report["frame"]  = frame;
    report["layer"]  = layer_str(layer);
    json targets_j = json::array();
    std::vector<std::string> warnings;
    bool all_regions_matched = true;

    for (ModelVolume *mv : targets) {
        const indexed_triangle_set &its = mv->mesh().its;
        TriangleSelector selector(mv->mesh());
        // Fresh canvas — regions declaratively define the final state; prior
        // paint on this volume (in the targeted layer) is discarded.
        selector.reset();

        // World-frame transform: predicates and reachable_from seed are in
        // world coords, so we transform centroids/normals BEFORE evaluating.
        // Written state still lands on mesh-local triangle indices (unchanged
        // by the coordinate reinterpretation).
        Transform3d trafo = Transform3d::Identity();
        if (frame_is_world) {
            trafo = mv->get_matrix();
            if (mv->get_object() && !mv->get_object()->instances.empty())
                trafo = mv->get_object()->instances.front()->get_matrix() * trafo;
        }
        const Matrix3d normal_trafo = trafo.linear().inverse().transpose();

        // Precompute per-triangle geometry once; regions reuse. Neighbors are
        // built lazily (only if a region actually needs flood-fill / CC).
        const size_t n_tris = its.indices.size();
        std::vector<Vec3d> centroids(n_tris);
        std::vector<Vec3d> normals(n_tris);
        std::vector<double> areas(n_tris);
        for (size_t i = 0; i < n_tris; ++i) {
            const auto &t = its.indices[i];
            Vec3d a = its.vertices[t(0)].cast<double>();
            Vec3d b = its.vertices[t(1)].cast<double>();
            Vec3d c = its.vertices[t(2)].cast<double>();
            if (frame_is_world) {
                a = trafo * a;
                b = trafo * b;
                c = trafo * c;
            }
            centroids[i] = (a + b + c) / 3.0;
            Vec3d cross  = (b - a).cross(c - a);
            areas[i]     = 0.5 * cross.norm();
            normals[i]   = areas[i] > 1e-12 ? (cross / (2.0 * areas[i])).eval() : Vec3d::UnitZ();
            // For scaled/mirrored volumes the cross-product normal is already
            // in world orientation; normal_trafo is only needed when we later
            // want a mesh-local normal — not the case here.
            (void)normal_trafo;
        }
        std::optional<FaceNeighbors> neighbors; // lazy
        std::optional<AABBTreeIndirect::Tree<3, float>> mesh_tree; // lazy (shape path)

        json volume_report;
        volume_report["object_name"]   = mv->get_object() ? mv->get_object()->name : "";
        volume_report["volume_name"]   = mv->name;
        volume_report["n_facets_total"] = n_tris;

        json regions_j = json::array();
        for (size_t ri = 0; ri < regions.size(); ++ri) {
            const Region &reg = regions[ri];

            // ── Shape-splat path (line/rect/polygon) ─────────────────
            // Projects a 2D shape onto the mesh via raycasting and splats
            // a Circle cursor at each hit, subdividing triangles at the
            // cursor edge. Bypasses the whole-triangle predicate below —
            // shape and whole-triangle predicates are mutually exclusive.
            if (reg.shape.kind != Shape::Kind::None) {
                if (frame_is_world) {
                    // Shape coords are always in the mesh-local frame for now;
                    // world-frame + shape is a follow-up (needs shape-side
                    // trafo application, easy but out of scope for MVP).
                    warnings.push_back("region[" + std::to_string(ri) +
                                       "] shape ignored: frame=world + shape not implemented");
                }
                if (!mesh_tree)
                    mesh_tree.emplace(AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
                        its.vertices, its.indices));

                const Shape &sh = reg.shape;
                const BoundingBoxf sbb = shape_bbox_2d(sh);

                // Ray direction along projection axis. Origin sits well
                // outside the mesh so the whole model is in front of the ray.
                const BoundingBoxf3 mesh_bb = bounding_box(its);
                Vec3d ray_dir;
                double ray_origin_z;
                if (sh.project == ProjectAxis::NegZ) {
                    ray_dir      = Vec3d(0.0, 0.0, -1.0);
                    ray_origin_z = mesh_bb.max.z() + 10.0;
                } else {
                    ray_dir      = Vec3d(0.0, 0.0, +1.0);
                    ray_origin_z = mesh_bb.min.z() - 10.0;
                }

                // Sample grid.
                const double  spacing = sh.sample_spacing_mm;
                const int     nx = std::max(1, int(std::ceil((sbb.max.x() - sbb.min.x()) / spacing)) + 1);
                const int     ny = std::max(1, int(std::ceil((sbb.max.y() - sbb.min.y()) / spacing)) + 1);

                // Fake "camera" for the cursor — placed at ray origin. Same
                // world coords as ray, since we're already in mesh space.
                const Transform3d identity = Transform3d::Identity();

                size_t samples_in_shape = 0, samples_hit = 0, samples_miss = 0, samples_clipped = 0;
                BoundingBoxf3 hit_bb;  bool hit_bb_have = false;

                for (int iy = 0; iy < ny; ++iy) {
                    const double y = sbb.min.y() + iy * spacing;
                    for (int ix = 0; ix < nx; ++ix) {
                        const double x = sbb.min.x() + ix * spacing;
                        const Vec2d p2(x, y);
                        if (!point_in_shape_2d(sh, p2)) continue;
                        ++samples_in_shape;

                        // Iterate hits along the ray until we find one whose
                        // point isn't on the clipped side. Without a clip
                        // the first hit wins. Guard the loop; a mesh with
                        // 8 stacked surfaces along one ray is pathological.
                        const Vec3d origin(x, y, ray_origin_z);
                        Vec3d origin_used = origin;
                        igl::Hit<float> hit;
                        bool hit_valid = AABBTreeIndirect::intersect_ray_first_hit(
                                its.vertices, its.indices, *mesh_tree, origin_used, ray_dir, hit);
                        Vec3d hit_pos_d;
                        if (hit_valid) hit_pos_d = origin_used + double(hit.t) * ray_dir;
                        if (hit_valid && clip.is_active()) {
                            int guard = 0;
                            while (hit_valid && clip.is_mesh_point_clipped(hit_pos_d.cast<float>()) && ++guard < 8) {
                                origin_used = hit_pos_d + 1e-3 * ray_dir;
                                hit_valid = AABBTreeIndirect::intersect_ray_first_hit(
                                        its.vertices, its.indices, *mesh_tree, origin_used, ray_dir, hit);
                                if (hit_valid) hit_pos_d = origin_used + double(hit.t) * ray_dir;
                            }
                            if (hit_valid && clip.is_mesh_point_clipped(hit_pos_d.cast<float>())) {
                                ++samples_clipped;
                                continue;
                            }
                        }
                        if (!hit_valid) {
                            ++samples_miss;
                            continue;
                        }
                        ++samples_hit;
                        const Vec3f hit_pos   = hit_pos_d.cast<float>();
                        const Vec3f origin_f  = origin_used.cast<float>();
                        if (!hit_bb_have) { hit_bb.min = hit_bb.max = hit_pos_d; hit_bb_have = true; }
                        else               hit_bb.merge(hit_pos_d);

                        auto cursor = TriangleSelector::SinglePointCursor::cursor_factory(
                            hit_pos, /*camera_pos*/ origin_f,
                            float(sh.cursor_radius_mm),
                            TriangleSelector::CursorType::CIRCLE,
                            identity, clip);
                        selector.select_patch(hit.id, std::move(cursor),
                                              reg.state, identity,
                                              /*triangle_splitting=*/true);
                    }
                }

                // Report — shape path uses samples/hits, not facets_matched;
                // provide both for consistency with the other regions.
                json region_j;
                region_j["index"]          = ri;
                region_j["kind"]           = reg.kind_str;
                region_j["shape"]          = (sh.kind == Shape::Kind::Line    ? "line"    :
                                              sh.kind == Shape::Kind::Rect    ? "rect"    :
                                              sh.kind == Shape::Kind::Polygon ? "polygon" : "?");
                region_j["samples_in_shape"] = samples_in_shape;
                region_j["samples_hit"]      = samples_hit;
                region_j["samples_miss"]     = samples_miss;
                if (samples_clipped > 0) region_j["samples_clipped"] = samples_clipped;
                region_j["facets_matched"]   = samples_hit;   // for --strict + summary parity
                if (hit_bb_have) region_j["bbox"] = bbox_to_json(hit_bb);

                if (samples_in_shape == 0) {
                    all_regions_matched = false;
                    warnings.push_back("region[" + std::to_string(ri) +
                                       "] shape: 0 samples inside — check sample_spacing_mm vs shape size");
                } else if (samples_hit == 0) {
                    all_regions_matched = false;
                    warnings.push_back("region[" + std::to_string(ri) +
                                       "] shape's projection fell entirely outside the mesh");
                }
                regions_j.push_back(std::move(region_j));
                continue;   // skip whole-triangle path
            }

            // Phase 1 — evaluate the field predicates over every source triangle.
            std::vector<bool> matched(n_tris, false);
            for (size_t i = 0; i < n_tris; ++i)
                if (predicate_matches(reg.pred, centroids[i], normals[i])) matched[i] = true;

            // Phase 2 — if reachable_from, flood-fill from the nearest triangle,
            // gated by the base-match predicate. Any base-matched triangle that
            // is NOT connected to the seed drops out.
            std::string fill_note;
            if (reg.pred.reachable_from) {
                if (!neighbors) neighbors.emplace(compute_face_neighbors(its));
                int seed = find_seed_triangle(its, *reg.pred.reachable_from);
                if (seed < 0 || !matched[seed]) {
                    fill_note = "reachable_from seed did not match the region's other predicates; nothing painted";
                    std::fill(matched.begin(), matched.end(), false);
                } else {
                    std::vector<bool> flood(n_tris, false);
                    std::queue<int> q;
                    q.push(seed); flood[seed] = true;
                    while (!q.empty()) {
                        int t = q.front(); q.pop();
                        for (int nb : (*neighbors)[t])
                            if (nb >= 0 && matched[nb] && !flood[nb]) { flood[nb] = true; q.push(nb); }
                    }
                    matched = std::move(flood);
                }
            }

            // Phase 3 — if connected_component, keep only the largest.
            if (reg.pred.connected_component && !reg.pred.reachable_from) {
                if (!neighbors) neighbors.emplace(compute_face_neighbors(its));
                matched = largest_connected_component(matched, *neighbors);
            }

            // Commit + collect stats.
            size_t        n_matched = 0;
            BoundingBoxf3 hit_bb;
            bool          hit_bb_have = false;
            double        area_sum    = 0.0;
            for (size_t i = 0; i < n_tris; ++i) {
                if (!matched[i]) continue;
                selector.set_facet(int(i), reg.state);
                ++n_matched;
                area_sum += areas[i];
                const auto &t = its.indices[i];
                for (int k = 0; k < 3; ++k) {
                    const Vec3d v = its.vertices[t(k)].cast<double>();
                    if (!hit_bb_have) { hit_bb.min = hit_bb.max = v; hit_bb_have = true; }
                    else               hit_bb.merge(v);
                }
            }

            json region_j;
            region_j["index"]          = ri;
            region_j["kind"]           = reg.kind_str;
            region_j["facets_matched"] = n_matched;
            region_j["area_mm2"]       = area_sum;
            if (hit_bb_have) region_j["bbox"] = bbox_to_json(hit_bb);
            if (!fill_note.empty()) region_j["note"] = fill_note;
            regions_j.push_back(std::move(region_j));

            if (n_matched == 0) {
                all_regions_matched = false;
                warnings.push_back("region[" + std::to_string(ri) + "] (" +
                                   reg.kind_str + ") matched 0 facets on volume \"" +
                                   mv->name + "\"");
            }
        }
        volume_report["regions"] = std::move(regions_j);

        // Commit into the targeted paint layer, replacing prior state there.
        FacetsAnnotation &fa = layer_facets(*mv, layer);
        fa.reset();
        fa.set(selector);

        targets_j.push_back(std::move(volume_report));
    }
    report["targets"]  = std::move(targets_j);
    report["warnings"] = warnings;

    report_out << report.dump(2) << std::endl;
    return all_regions_matched;
}

// ================================ --render-paint ============================

namespace {

// Small headless PPM writer. P6 (binary), 3-byte RGB. Standard image
// viewers and ImageMagick accept it; convert with `magick out.ppm out.png`.
bool write_ppm(const std::string &path, int w, int h, const std::vector<uint8_t> &rgb)
{
    boost::nowide::ofstream f(path.c_str(), std::ios::binary);
    if (!f) return false;
    f << "P6\n" << w << " " << h << "\n255\n";
    f.write(reinterpret_cast<const char*>(rgb.data()), std::streamsize(rgb.size()));
    return f.good();
}

// Color palette for MMU extruders (1–16). ENFORCER=red, BLOCKER=blue,
// FUZZY_SKIN=green, unpainted=lambert gray. Miss=black.
std::array<uint8_t, 3> color_for_state(EnforcerBlockerType st, Layer layer)
{
    if (layer == Layer::Supports || layer == Layer::Seam) {
        if (st == EnforcerBlockerType::ENFORCER) return {240,  60,  40};   // red
        if (st == EnforcerBlockerType::BLOCKER)  return { 40,  80, 240};   // blue
    }
    if (layer == Layer::Fuzzy && st == EnforcerBlockerType::FUZZY_SKIN)
        return { 60, 200,  80};   // green
    if (layer == Layer::Mmu) {
        static constexpr std::array<std::array<uint8_t,3>, 16> palette = {{
            {230, 60, 60}, {60,120,230}, {60,190, 90}, {230,180, 60},
            {200, 60,180}, {60,200,220}, {230,110, 40}, {150, 90,210},
            {200,220, 60}, {120,220,180},{220, 60,120},{110,180, 40},
            { 90,110,200}, {230,140,110},{190,100,140},{180,180,190}
        }};
        int idx = int(st) - 1;
        if (idx >= 0 && idx < 16) return palette[idx];
    }
    return {180, 180, 180};   // unknown / default
}

} // namespace

bool render_paint(const Model &model, const std::string &out_path,
                  const RenderPaintOpts &opts, std::ostream &report_out)
{
    Layer layer;
    try { layer = parse_layer(opts.layer); }
    catch (const std::exception &e) {
        report_out << json({{ "error", "invalid_layer" }, { "detail", e.what() }}).dump(2) << std::endl;
        return false;
    }
    if (opts.view != "top" && opts.view != "bottom") {
        report_out << json({{ "error", "invalid_view" }, { "detail", "view must be top|bottom" }}).dump(2) << std::endl;
        return false;
    }
    const bool from_top = (opts.view == "top");

    TriangleSelector::ClippingPlane clip{};
    if (std::isfinite(opts.clip_z_below)) {
        clip.normal = { 0.f, 0.f, 1.f };
        clip.offset = float(opts.clip_z_below);
    } else if (std::isfinite(opts.clip_z_above)) {
        clip.normal = { 0.f, 0.f, -1.f };
        clip.offset = float(-opts.clip_z_above);
    }

    // World bbox over all objects × instance 0. Sets the camera frame.
    BoundingBoxf3 world_bb;
    bool have_bb = false;
    for (const ModelObject *mo : model.objects) {
        if (!mo || mo->instances.empty()) continue;
        for (size_t vi = 0; vi < mo->volumes.size(); ++vi) {
            const ModelVolume *mv = mo->volumes[vi];
            if (!mv || !mv->is_model_part()) continue;
            const indexed_triangle_set &its = mv->mesh().its;
            Transform3d trafo = mo->instances.front()->get_matrix() * mv->get_matrix();
            for (const Vec3f &v : its.vertices) {
                Vec3d w = trafo * v.cast<double>();
                if (!have_bb) { world_bb.min = world_bb.max = w; have_bb = true; }
                else            world_bb.merge(w);
            }
        }
    }
    if (!have_bb) {
        report_out << json({{ "error", "no_geometry" }}).dump(2) << std::endl;
        return false;
    }

    // Pad the frame slightly so the model isn't flush with the edges.
    Vec3d size = world_bb.max - world_bb.min;
    Vec3d pad  = 0.05 * size;
    world_bb.min -= pad;
    world_bb.max += pad;

    const int W = std::max(64, opts.width_px);
    const int H = std::max(64, opts.height_px);
    std::vector<uint8_t> rgb(size_t(W) * size_t(H) * 3, 0);

    // Aspect-preserved pixel → world XY.
    const double frame_w = world_bb.max.x() - world_bb.min.x();
    const double frame_h = world_bb.max.y() - world_bb.min.y();
    const double px_size = std::max(frame_w / W, frame_h / H);   // world mm per pixel

    // Ray direction & origin z.
    const Vec3d ray_dir = from_top ? Vec3d(0, 0, -1) : Vec3d(0, 0, 1);
    const double origin_z = from_top ? (world_bb.max.z() + 10.0) : (world_bb.min.z() - 10.0);

    // Precompute per (object, volume) the world-transformed itset + AABB tree
    // for the base mesh plus one itset+tree per painted state. Rays are cast
    // against all of these and the closest hit determines the pixel color.
    struct RenderVolume {
        indexed_triangle_set                     base_world;   // vertices transformed to world
        AABBTreeIndirect::Tree<3, float>         base_tree;
        std::vector<Vec3d>                       base_face_normals;
        struct StatePaint {
            EnforcerBlockerType state;
            indexed_triangle_set its;                          // world coords
            AABBTreeIndirect::Tree<3, float> tree;
        };
        std::vector<StatePaint> painted;
    };

    std::vector<RenderVolume> render_vols;
    for (const ModelObject *mo : model.objects) {
        if (!mo || mo->instances.empty()) continue;
        Transform3d obj_trafo = mo->instances.front()->get_matrix();
        for (const ModelVolume *mv : mo->volumes) {
            if (!mv || !mv->is_model_part()) continue;
            Transform3d trafo = obj_trafo * mv->get_matrix();

            RenderVolume rv;
            const indexed_triangle_set &src = mv->mesh().its;
            rv.base_world.indices  = src.indices;
            rv.base_world.vertices.reserve(src.vertices.size());
            for (const Vec3f &v : src.vertices)
                rv.base_world.vertices.emplace_back((trafo * v.cast<double>()).cast<float>());
            rv.base_tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
                rv.base_world.vertices, rv.base_world.indices);
            rv.base_face_normals.reserve(rv.base_world.indices.size());
            for (const auto &t : rv.base_world.indices) {
                const Vec3d a = rv.base_world.vertices[t(0)].cast<double>();
                const Vec3d b = rv.base_world.vertices[t(1)].cast<double>();
                const Vec3d c = rv.base_world.vertices[t(2)].cast<double>();
                Vec3d n = (b - a).cross(c - a);
                double norm = n.norm();
                Vec3d unit_n;
                if (norm > 1e-12) unit_n = n / norm;
                else              unit_n = Vec3d(0, 0, 1);
                rv.base_face_normals.emplace_back(unit_n);
            }

            // Painted itsets per state — from FacetsAnnotation. Coords come
            // back in mesh-local so we transform them here.
            const FacetsAnnotation &fa = const_cast<ModelVolume*>(mv)->supported_facets; // upstream lacks const overload
            const FacetsAnnotation *layer_fa = nullptr;
            switch (layer) {
            case Layer::Supports: layer_fa = &mv->supported_facets;         break;
            case Layer::Seam:     layer_fa = &mv->seam_facets;              break;
            case Layer::Mmu:      layer_fa = &mv->mmu_segmentation_facets;  break;
            case Layer::Fuzzy:    layer_fa = &mv->fuzzy_skin_facets;        break;
            }
            (void)fa;
            if (layer_fa && !layer_fa->empty()) {
                std::vector<std::pair<EnforcerBlockerType, std::string>> states =
                    (layer == Layer::Mmu)   ? mmu_states() :
                    (layer == Layer::Fuzzy) ? fuzzy_states() :
                                              supports_states();
                for (const auto &st : states) {
                    if (!layer_fa->has_facets(*mv, st.first)) continue;
                    indexed_triangle_set painted = layer_fa->get_facets_strict(*mv, st.first);
                    if (painted.indices.empty()) continue;
                    RenderVolume::StatePaint sp;
                    sp.state = st.first;
                    sp.its.indices = painted.indices;
                    sp.its.vertices.reserve(painted.vertices.size());
                    for (const Vec3f &v : painted.vertices)
                        sp.its.vertices.emplace_back((trafo * v.cast<double>()).cast<float>());
                    sp.tree = AABBTreeIndirect::build_aabb_tree_over_indexed_triangle_set(
                        sp.its.vertices, sp.its.indices);
                    rv.painted.push_back(std::move(sp));
                }
            }
            render_vols.push_back(std::move(rv));
        }
    }

    // Rasterize.
    for (int py = 0; py < H; ++py) {
        const double y = world_bb.max.y() - (py + 0.5) * px_size;
        for (int px = 0; px < W; ++px) {
            const double x = world_bb.min.x() + (px + 0.5) * px_size;
            const Vec3d origin(x, y, origin_z);

            // Track the closest hit across all (mesh, state itsets, volumes).
            // Painted-state hits are preferred when their t is within a small
            // epsilon of the base hit (they sit on the same plane).
            double best_t = std::numeric_limits<double>::infinity();
            std::array<uint8_t, 3> pixel = {0, 0, 0};    // black = miss

            auto try_hit = [&](const RenderVolume &rv, const AABBTreeIndirect::Tree<3, float> &tree,
                               const indexed_triangle_set &target, bool is_painted,
                               EnforcerBlockerType st) {
                // Advance past clipped hits (same idea as apply_spec).
                Vec3d origin_used = origin;
                igl::Hit<float> hit;
                bool ok = AABBTreeIndirect::intersect_ray_first_hit(
                            target.vertices, target.indices, tree, origin_used, ray_dir, hit);
                if (ok && clip.is_active()) {
                    int guard = 0;
                    Vec3d hp = origin_used + double(hit.t) * ray_dir;
                    while (ok && clip.is_mesh_point_clipped(hp.cast<float>()) && ++guard < 8) {
                        origin_used = hp + 1e-3 * ray_dir;
                        ok = AABBTreeIndirect::intersect_ray_first_hit(
                                target.vertices, target.indices, tree, origin_used, ray_dir, hit);
                        if (ok) hp = origin_used + double(hit.t) * ray_dir;
                    }
                    if (ok && clip.is_mesh_point_clipped(hp.cast<float>())) return;
                }
                if (!ok) return;
                // Total distance from ORIGINAL origin along ray_dir.
                const Vec3d hp = origin_used + double(hit.t) * ray_dir;
                const double t_world = (hp - origin).dot(ray_dir);

                if (is_painted) {
                    // Painted sub-tris sit ON the base tri plane → prefer over
                    // a base-only hit at the same depth.
                    if (t_world <= best_t + 1e-3) {
                        best_t = std::min(best_t, t_world);
                        pixel = color_for_state(st, layer);
                    }
                    return;
                }
                if (t_world < best_t) {
                    best_t = t_world;
                    // Lambert shading vs the -view direction so top-lit.
                    const Vec3d n = rv.base_face_normals[size_t(hit.id)];
                    double lit = std::max(0.15, -n.dot(ray_dir));
                    const uint8_t g = uint8_t(std::round(lit * 210.0));
                    pixel = { g, g, g };
                }
            };

            for (const RenderVolume &rv : render_vols) {
                try_hit(rv, rv.base_tree, rv.base_world, /*is_painted=*/false, EnforcerBlockerType::NONE);
                for (const auto &sp : rv.painted)
                    try_hit(rv, sp.tree, sp.its, /*is_painted=*/true, sp.state);
            }

            size_t off = (size_t(py) * W + size_t(px)) * 3;
            rgb[off + 0] = pixel[0];
            rgb[off + 1] = pixel[1];
            rgb[off + 2] = pixel[2];
        }
    }

    if (!write_ppm(out_path, W, H, rgb)) {
        report_out << json({{ "error", "cannot_write_ppm" }, { "path", out_path }}).dump(2) << std::endl;
        return false;
    }
    json report;
    report["out_path"]  = out_path;
    report["layer"]     = opts.layer;
    report["view"]      = opts.view;
    report["width_px"]  = W;
    report["height_px"] = H;
    report["px_size_mm"] = px_size;
    report["clip"] = json::object();
    if (std::isfinite(opts.clip_z_below)) report["clip"]["z_below"] = opts.clip_z_below;
    if (std::isfinite(opts.clip_z_above)) report["clip"]["z_above"] = opts.clip_z_above;
    report["frame_world_bbox"] = bbox_to_json(world_bb);
    report_out << report.dump(2) << std::endl;
    return true;
}

} // namespace PaintCLI
} // namespace Slic3r
