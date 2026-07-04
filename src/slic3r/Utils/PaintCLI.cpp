// PaintCLI.cpp — CLI paint-inspection primitives. See PaintCLI.hpp.
#include "PaintCLI.hpp"

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
                     "Paint gizmos and the future --paint-supports predicates "
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

// =============================== --paint-supports ===========================

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

struct Region {
    std::string         kind_str; // for report
    EnforcerBlockerType state;    // ENFORCER / BLOCKER / NONE (=clear)
    Predicate           pred;
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

EnforcerBlockerType parse_kind(const std::string &s)
{
    if (s == "enforcer") return EnforcerBlockerType::ENFORCER;
    if (s == "blocker")  return EnforcerBlockerType::BLOCKER;
    if (s == "clear")    return EnforcerBlockerType::NONE;
    throw std::runtime_error("--paint-supports: unknown region kind \"" + s +
                             "\" (must be enforcer|blocker|clear)");
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

bool apply_supports_spec(Model &model, const std::string &spec_path,
                         const std::string &source_path,
                         std::ostream &report_out)
{
    json spec;
    {
        boost::nowide::ifstream in(spec_path.c_str());
        if (!in) {
            BOOST_LOG_TRIVIAL(error) << "--paint-supports: cannot open " << spec_path;
            report_out << json({{ "error", "cannot_open_spec" }, { "path", spec_path }}).dump(2) << std::endl;
            return false;
        }
        try { in >> spec; }
        catch (const std::exception &e) {
            BOOST_LOG_TRIVIAL(error) << "--paint-supports: invalid JSON in " << spec_path << ": " << e.what();
            report_out << json({{ "error", "invalid_json" }, { "path", spec_path }, { "detail", e.what() }}).dump(2) << std::endl;
            return false;
        }
    }

    // Parse regions upfront so errors surface before we touch the model.
    std::vector<Region> regions;
    if (!spec.contains("regions") || !spec["regions"].is_array()) {
        report_out << json({{ "error", "regions_missing_or_not_array" }}).dump(2) << std::endl;
        return false;
    }
    for (const json &r : spec["regions"]) {
        Region reg;
        reg.kind_str = r.value("kind", "");
        reg.state    = parse_kind(reg.kind_str);
        reg.pred     = parse_predicate(r.value("select", json::object()));
        regions.push_back(std::move(reg));
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
    report["frame"]  = "mesh_local";
    json targets_j = json::array();
    std::vector<std::string> warnings;
    bool all_regions_matched = true;

    for (ModelVolume *mv : targets) {
        const indexed_triangle_set &its = mv->mesh().its;
        TriangleSelector selector(mv->mesh());
        // Fresh canvas — regions declaratively define the final state; prior
        // paint on this volume is discarded (matches "spec-driven" semantics).
        selector.reset();

        // Precompute per-triangle geometry once; regions reuse. Neighbors are
        // built lazily (only if a region actually needs flood-fill / CC).
        const size_t n_tris = its.indices.size();
        std::vector<Vec3d> centroids(n_tris);
        std::vector<Vec3d> normals(n_tris);
        std::vector<double> areas(n_tris);
        for (size_t i = 0; i < n_tris; ++i) {
            const auto &t = its.indices[i];
            const Vec3d a = its.vertices[t(0)].cast<double>();
            const Vec3d b = its.vertices[t(1)].cast<double>();
            const Vec3d c = its.vertices[t(2)].cast<double>();
            centroids[i] = (a + b + c) / 3.0;
            Vec3d cross  = (b - a).cross(c - a);
            areas[i]     = 0.5 * cross.norm();
            normals[i]   = areas[i] > 1e-12 ? (cross / (2.0 * areas[i])).eval() : Vec3d::UnitZ();
        }
        std::optional<FaceNeighbors> neighbors; // lazy

        json volume_report;
        volume_report["object_name"]   = mv->get_object() ? mv->get_object()->name : "";
        volume_report["volume_name"]   = mv->name;
        volume_report["n_facets_total"] = n_tris;

        json regions_j = json::array();
        for (size_t ri = 0; ri < regions.size(); ++ri) {
            const Region &reg = regions[ri];

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

        // Commit selector state into the FacetsAnnotation, replacing prior paint.
        mv->supported_facets.reset();
        mv->supported_facets.set(selector);

        targets_j.push_back(std::move(volume_report));
    }
    report["targets"]  = std::move(targets_j);
    report["warnings"] = warnings;

    report_out << report.dump(2) << std::endl;
    return all_regions_matched;
}

} // namespace PaintCLI
} // namespace Slic3r
