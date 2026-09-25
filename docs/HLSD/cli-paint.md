# CLI Painting: High Level Design

## Why it exists

The paint gizmos — supports, seam, MMU colour, fuzzy skin — store per-facet state on the
model, and until now the only way to set that state was to open the GUI and drag a brush
over the part. A scripted or automated pipeline that wanted a support enforcer on one
overhang had to hand the job back to a person.

`--paint` applies a JSON spec of ordered regions to one paint layer of the loaded model,
in place. Chained with `--export-3mf` it saves a painted project; chained with `--slice`
it bakes the paint into G-code in the same run. `--inspect-paint` is the read side: it
reports what is already painted, so a pipeline can inspect, decide and write without a
display.

## Two ways a region selects facets

A region's `select` object is read by both selector families, and which one runs depends
on the keys present.

**Whole source triangles.** `bbox`, `z_min`, `z_max` and `normal_below` match triangles of
the mesh as it was loaded; `reachable_from` flood-fills from the triangle nearest a point,
and `connected_component` keeps only the largest connected subset of the match. Every
matched triangle is painted entirely. On a coarse mesh this is all-or-nothing: a box lid
that is two triangles is either wholly painted or wholly missed.

**A projected 2D shape.** `line {a, b, width_mm}`, `rect {min, max}` and
`polygon {points}` describe a figure in the projection plane. The figure is sampled on a
grid of `sample_spacing_mm`, each sample is cast along `project` (`-Z` by default, `+Z`
also accepted), and where a ray hits the mesh the hit is splatted with a circular cursor
of `cursor_radius_mm` through `TriangleSelector::select_patch()` with triangle splitting
enabled.

That last part is the point of the shape form: the selector subdivides the triangles under
the cursor, exactly as the GUI gizmos do, so a shape can paint a small patch of a large
facet. A 4 x 4 mm rectangle on a part whose target facet is 414 mm² paints under 4% of
that one triangle, turning 256 source facets into ~11500 sub-facets in the stored
annotation. The whole-triangle selectors cannot express that, and no combination of them
approximates it.

## Frames

`frame` selects the coordinate system the spec is written in: `mesh_local` (the default)
or `world`. Mesh-local is the *centred* frame — OrcaSlicer centres a mesh when it loads
it — and it is the frame `--inspect-mesh` prints as `bbox_mesh_local`, so those numbers
can be used directly in a spec.

Shapes are implemented for `mesh_local` only. A shape given with `frame: "world"` is
skipped with the warning `shape ignored: frame=world + shape not implemented`, and the
region matches nothing. The whole-triangle selectors honour both frames.

## Reporting

Every run prints a JSON coverage report to stdout: per region, the shape kind, how many
samples fell inside the figure, how many hit the mesh, how many missed, and the bounding
box of the hits. The report is the contract a caller checks; under `--strict` a region
that matched zero facets fails the action rather than leaving the caller to notice an
empty result.

## Constraints worth knowing

A shape that reports many `samples_in_shape` but `samples_hit: 0`, with
`shape's projection fell entirely outside the mesh`, usually sits over a **hole** in the
part rather than off the model. The message names the projection, not the coordinates, and
a part that is hollow in the middle produces it for a figure centred on the origin while a
figure a few millimetres to the side hits every sample. Pick the target from the geometry —
the centroid of a large upward-facing facet, say — rather than from the bounding box.

`cursor_radius_mm` spills past the edge of the figure, and where the figure reaches the
edge of a face the cursor wraps around onto the adjacent one. A rectangle painted at the
rim of a part therefore reports a bounding box taller than the figure.

`kind` is `enforcer`, `blocker` or `clear` for supports and seam, `fuzzy_skin` for the
fuzzy layer, and `extruder_N` for MMU.
