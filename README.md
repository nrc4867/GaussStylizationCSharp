# GaussStylization (C#)

A managed C# library for **Gauss Stylization**: it deforms a triangle mesh so that its face normals
gather around a set of *preferred directions* while the surface stays as rigid as possible.

This is a C# port of [ugogon/gaussStylization](https://github.com/ugogon/gaussStylization) and is
licensed under the **MPL-2.0** like the original. See [`NOTICE.md`](NOTICE.md) for attribution and
[`LICENSE`](LICENSE) for the license text.

* Pure managed code (`netstandard2.1`), no native dependencies and no other packages.
* Plain API: arrays in, arrays out.
* No UI; just the algorithm and its tests.

## Usage

```csharp
using GaussStylization;

double[] vertices  = ...; // x0,y0,z0, x1,y1,z1, ...   (3 * vertexCount)
int[]    triangles = ...; // a0,b0,c0, a1,b1,c1, ...   (3 * triangleCount)

double[] stylized = GaussStylizer.Stylize(
    vertices,
    triangles,
    preferredDirections: GaussStylizer.CubeDirections(), // xyz triples, any number of them
    vertexWeights: null,                                 // optional, one weight per vertex
    gridAxes: null,                                      // optional 3x3 frame, row-major
    iterations: 50);                                     // outer iterations

// stylized has the same layout as vertices; the inputs are not modified.
```

| Parameter | Meaning |
|---|---|
| `vertices` | Vertex positions, flat xyz. |
| `triangles` | Triangle vertex indices, flat. The mesh must be edge-manifold; open meshes and several components are fine. |
| `preferredDirections` | The preferred surface normals, flat xyz (normalised internally, must be distinct). Expressed in the grid frame. |
| `vertexWeights` | Optional, per vertex, usually in `[0,1]`. A face is stylized with the mean weight of its vertices; `0` leaves a region alone. `null` = 1 everywhere. |
| `gridAxes` | Optional orientation of the grid: three orthonormal axis vectors `[ax,ay,az, bx,by,bz, cx,cy,cz]`. A preferred direction `(d0,d1,d2)` points along `d0*A + d1*B + d2*C` in world space. `null` = world axes. Rotating the axes rotates all preferred directions (the rotate keys of the original app). |
| `iterations` | Outer iterations; each is one local step plus one global solve. `0` returns a copy of the input. |
| `options` | Optional `GaussStylizationOptions`: `Lambda` (4), `InnerIterations` (1), `FixedVertices`, solver tolerance. Defaults match the original application. |

A `float[]` overload exists as well; it computes in double precision.

### Styles: rings and style groups

The call above uses one style made of discrete directions. For cylinders and cones, per-style parameters,
or different styles on different parts of the mesh, pass `GaussStyle` objects and an optional style index per triangle:

```csharp
var styles = new[]
{
    GaussStyle.Cube(),                                  // group 0
    new GaussStyle(new double[] { 0, 1, 0, 0, -1, 0 })  // group 1: caps and a ring = cylinder
    {
        RingHeights = new[] { 0.0 },
        Mu = 1, Sigma = 4, RingContribution = 0.5,
    },
};
int[] faceGroups = ...; // one entry per triangle, index into styles (null = all 0)

double[] stylized = GaussStylizer.Stylize(vertices, triangles, styles, faceGroups,
    vertexWeights: null, gridAxes: null, iterations: 50);
```

| `GaussStyle` member | Meaning |
|---|---|
| `Directions` | Discrete preferred normals (flat xyz), in the grid frame. May be empty if rings are given. |
| `RingHeights` | Optional rings: each height `d` in `[-1,1]` allows every normal with `n.y = d` (local Y of the grid). A ring at `0` allows all normals perpendicular to the axis. |
| `Mu`, `Sigma` | Strength and sharpness of this style (defaults 1 and 4). |
| `RingContribution` | Weight of the discrete directions relative to the rings (default 0.5; only used with rings). |
| `Cube()`, `Cylinder()`, `Cone()` | The presets of the original application. |

The ring axis is the grid's second axis, i.e. world Y unless `gridAxes` is given.

### Notes

* One vertex per connected component is held in place (otherwise the global solve is singular).
  Pass `options.FixedVertices` to choose which vertices stay put.
* Stylization is deterministic. Meshes with at least 1000 faces/edges/vertices use multiple threads.
* Styles whose targets are rank deficient (for example a single ring, or only +-Y) make the per-face Newton step
  ambiguous, exactly as in the original; results for such styles are valid but can differ in tiny numerical
  details between rotated copies of the same problem.

## Build and test

```
dotnet build
dotnet test
dotnet pack GaussStylizationCSharp -c Release -o artifacts   # creates GaussStylization.0.1.0.nupkg locally
```

Nothing is published by this repository; `dotnet pack` only produces a local package.

## License

Mozilla Public License 2.0 — see [`LICENSE`](LICENSE) and [`NOTICE.md`](NOTICE.md).
