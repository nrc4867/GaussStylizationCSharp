# Notice and attribution

**GaussStylization (C#)** is a derivative work of
**[gaussStylization](https://github.com/ugogon/gaussStylization)** by Ugo Finnendahl
(commit `d73e129fc3d3faa53f17f53078093a860a5c45b3`), the reference C++ implementation of
*Gauss Stylization*. The original project is licensed under the
[Mozilla Public License 2.0](https://mozilla.org/MPL/2.0/); this port is distributed under the
same license (see [`LICENSE`](LICENSE)). The full license text of the original project
(`LICENSE.MPL2`) is reproduced unmodified in `LICENSE`.

## What is derived from what

| This library | Derived from |
|---|---|
| `Internal/GaussStylizationCore.cs` (face / edge / local / global steps) | `utils_gauss/gauss_style_single_iteration.cpp`, `utils_gauss/gauss_style_precomputation.cpp` (gaussStylization) |
| `Internal/GFunction.cs` (g function incl. rings, gradient, Hessian, weight normalisation) | `utils_gauss/g_function.cpp`, `utils_gauss/normalize_g.cpp` (gaussStylization) |
| Rotation fit (`Mat3.FitRotation`) | `utils_cubic/orthogonal_procrustes.cpp` (Cubic Stylization code, bundled in gaussStylization) |
| Edge topology, cotangent weights, ARAP right-hand side (`MeshTopology`, `GaussStylizationCore.ArapRhs`) | [libigl](https://github.com/libigl/libigl): `edge_topology`, `cotmatrix_entries`, `arap_rhs` / `arap_linear_block`, `per_face_normals` (MPL-2.0, Copyright Alec Jacobson and others) |
| Fixed-vertex Laplacian solve (`LaplacianSolver`) | Replaces libigl `min_quad_with_fixed` / Eigen sparse factorisation with a managed conjugate-gradient solver |

The gaussStylization repository bundles code from *Cubic Stylization* by Hsueh-Ti Derek Liu and
Alec Jacobson, and builds on libigl. The credit for the algorithms belongs to those authors.

## Differences from the original

* Pure managed C#, no native dependencies (libigl, Eigen and ImGui are not used).
* The ring axis is not hard-coded to world Y: it is the second axis of the optional grid frame (world Y by default).
* Style groups are given as a per-triangle index into a list of `GaussStyle` objects instead of being painted in a GUI.
* A per-vertex weight scales how strongly each face is stylized (mean of the face's three vertex weights);
  the original selects discrete style groups per face instead.
* One vertex per connected component is held fixed (the original fixes a single vertex).
* The linear systems are solved iteratively (Jacobi-preconditioned CG) instead of by sparse Cholesky/LDLT.
* No normalisation of the input mesh to a unit box is performed; the cotangent formulation is scale invariant.

## Original work

* Gauss Stylization — <https://github.com/ugogon/gaussStylization> (MPL-2.0)
* Cubic Stylization — <https://www.dgp.toronto.edu/projects/cubic-stylization/>
* libigl — <https://libigl.github.io/> (MPL-2.0)

The example meshes of the original repository (CC BY, see its `license_meshes.txt`) are **not** part of
this project.

Every source file in this repository carries the MPL-2.0 header and refers to this notice.
