// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.
// Port of utils_gauss/gauss_style_precomputation.cpp and
// utils_gauss/gauss_style_single_iteration.cpp (single style group, discrete g only).

using System;
using System.Collections.Concurrent;
using System.Threading.Tasks;

namespace GaussStylization.Internal
{
    /// <summary>
    /// The Gauss stylization optimisation: alternating local steps (per-vertex rotation fit, per-face
    /// normal update and per-edge target update via an augmented Lagrangian) and a global Poisson-like
    /// solve with the cotangent Laplacian.
    /// </summary>
    internal sealed class GaussStylizationCore
    {
        private readonly double[] _v;            // original positions (3n)
        private readonly MeshTopology _topo;
        private readonly GFunction[] _g;     // per style group
        private readonly double[] _groupMu;     // per style group
        private readonly int[] _faceGroup;
        private readonly double[] _faceStrength; // per-face multiplier on g and mu
        private readonly double _lambda;
        private readonly int _inner;
        private readonly bool[] _fixed;
        private readonly LaplacianSolver _solver;

        public GaussStylizationCore(
            double[] vertices, int[] faces, GFunction[] groupFunctions, double[] groupMu, int[] faceGroups,
            double[]? vertexWeights, GaussStylizationOptions options, bool[] fixedVertices)
        {
            _v = vertices;
            _topo = new MeshTopology(vertices, faces);
            _g = groupFunctions;
            _groupMu = groupMu;
            _faceGroup = faceGroups;
            _lambda = options.Lambda;
            _inner = options.InnerIterations;
            _fixed = fixedVertices;
            _solver = new LaplacianSolver(_topo, fixedVertices, options.SolverTolerance, options.SolverMaxIterations);

            _faceStrength = new double[_topo.FaceCount];
            for (int f = 0; f < _topo.FaceCount; f++)
            {
                if (vertexWeights == null) { _faceStrength[f] = 1.0; continue; }
                _faceStrength[f] = (vertexWeights[faces[3 * f]] + vertexWeights[faces[3 * f + 1]] + vertexWeights[faces[3 * f + 2]]) / 3.0;
            }
        }

        public double[] Run(int iterations)
        {
            var u = (double[])_v.Clone();
            for (int i = 0; i < iterations; i++) SingleIteration(u);
            return u;
        }

        private void SingleIteration(double[] u)
        {
            int nf = _topo.FaceCount, ne = _topo.EdgeCount, nv = _topo.VertexCount;
            int[] ev = _topo.EdgeVertices, ef = _topo.EdgeFaces, fe = _topo.FaceEdges;

            Vec3[] nStar = _topo.FaceNormals(u);
            var eStar = new Vec3[ne];
            var uf = new double[3 * nf];

            ParallelFor(ne, e =>
            {
                eStar[e] = MeshTopology.Position(u, ev[2 * e + 1]) - MeshTopology.Position(u, ev[2 * e]);
            });

            for (int iter = 0; iter < _inner; iter++)
            {
                // ---- face step: update normals by a (projected) Newton step ------------------
                ParallelFor(nf, f =>
                {
                    double s = _faceStrength[f];
                    Vec3 n = nStar[f];

                    // gradient / Hessian of g w.r.t. n (equations 8, 9), scaled by the face weight
                    _g[_faceGroup[f]].GradientAndHessian(n, out Vec3 gg, out Mat3 hg);
                    gg = s * gg;
                    hg = s * hg;

                    Vec3 gr = new Vec3(0, 0, 0);
                    Mat3 hr = default;
                    for (int j = 0; j < 3; j++)
                    {
                        int e = fe[3 * f + j];
                        double w = _topo.EdgeWeight[e];
                        Vec3 es = eStar[e];
                        hr = hr + (_lambda * w) * Mat3.Outer(es, es);
                        gr = gr + (_lambda * w * (n.Dot(es) + uf[3 * f + j])) * es;
                    }

                    Vec3 gn = gr - gg;
                    Vec3 newton = Mat3.SolveSymmetric(hr - hg, -gn);

                    // project out the component along n
                    Vec3 d = newton - n.Dot(newton) * n;
                    if (d.Dot(gn) > 0)
                    {
                        // Newton points uphill: fall back to a damped gradient step
                        d = -0.1 * (gn - n.Dot(gn) * n);
                    }

                    nStar[f] = (n + d).Normalized();
                });

                // ---- edge step: closed-form update of the edge targets ------------------------
                ParallelFor(ne, ij =>
                {
                    Vec3 eij = MeshTopology.Position(u, ev[2 * ij + 1]) - MeshTopology.Position(u, ev[2 * ij]);
                    Mat3 a = Mat3.Identity;
                    Vec3 rhs = eij;
                    for (int side = 0; side < 2; side++)
                    {
                        int f = ef[2 * ij + side];
                        if (f == -1) continue;
                        double mu = _groupMu[_faceGroup[f]] * _faceStrength[f];
                        Vec3 n = nStar[f];
                        a = a + mu * Mat3.Outer(n, n);
                        for (int k = 0; k < 3; k++)
                            if (fe[3 * f + k] == ij) rhs = rhs - (mu * uf[3 * f + k]) * n;
                    }
                    eStar[ij] = Mat3.SolveSymmetric(a, rhs);
                });

                // ---- dual update --------------------------------------------------------------
                ParallelFor(nf, f =>
                {
                    for (int j = 0; j < 3; j++) uf[3 * f + j] += eStar[fe[3 * f + j]].Dot(nStar[f]);
                });
            }

            // right-hand side contribution of the target edges (equation 16)
            var target = new Vec3[nv];
            for (int e = 0; e < ne; e++)
            {
                double w = _topo.EdgeWeight[e];
                int v1 = ev[2 * e], v2 = ev[2 * e + 1];
                target[v1] = target[v1] - w * eStar[e];
                target[v2] = target[v2] + w * eStar[e];
            }

            // ---- local step: best-fit rotation per vertex -----------------------------------
            var rot = new Mat3[nv];
            ParallelFor(nv, i =>
            {
                Mat3 sb = default;
                for (int k = _topo.VertexFaceStart[i]; k < _topo.VertexFaceStart[i + 1]; k++)
                {
                    int f = _topo.VertexFaces[k];
                    for (int j = 0; j < 3; j++)
                    {
                        int e = fe[3 * f + j];
                        int a = ev[2 * e], b = ev[2 * e + 1];
                        Vec3 dv = MeshTopology.Position(_v, b) - MeshTopology.Position(_v, a);
                        Vec3 du = MeshTopology.Position(u, b) - MeshTopology.Position(u, a);
                        sb = sb + _topo.EdgeWeight[e] * Mat3.Outer(dv, du);
                    }
                }
                rot[i] = Mat3.FitRotation(sb);
            });

            // ---- global step -------------------------------------------------------------------
            var arap = ArapRhs(rot);
            var rhs = new double[3][];
            var x = new double[3][];
            var fixedValues = new double[3][];
            for (int d = 0; d < 3; d++)
            {
                rhs[d] = new double[nv];
                x[d] = new double[nv];
                fixedValues[d] = new double[nv];
            }
            for (int i = 0; i < nv; i++)
            {
                Vec3 bc = (arap[i] + _lambda * target[i]) * (1.0 / (1.0 + _lambda));
                rhs[0][i] = bc.X; rhs[1][i] = bc.Y; rhs[2][i] = bc.Z;
                for (int d = 0; d < 3; d++) { x[d][i] = u[3 * i + d]; fixedValues[d][i] = _v[3 * i + d]; }
            }
            Parallel.For(0, 3, d => _solver.Solve(rhs[d], fixedValues[d], x[d]));
            for (int i = 0; i < nv; i++)
                for (int d = 0; d < 3; d++) u[3 * i + d] = x[d][i];
        }

        /// <summary>
        /// Right-hand side B = K * R of the ARAP global step for the "spokes and rims" energy
        /// (igl::arap_rhs): per triangle edge (s,d) opposite vertex o, with p = cot/2 (V_s - V_d) / 3,
        /// B_s += (R_s + R_d + R_o) p and B_d -= (R_s + R_d + R_o) p.
        /// </summary>
        private Vec3[] ArapRhs(Mat3[] rot)
        {
            int nv = _topo.VertexCount;
            var b = new Vec3[nv];
            int[] faces = _topo.Faces;
            for (int f = 0; f < _topo.FaceCount; f++)
            {
                for (int j = 0; j < 3; j++)
                {
                    int o = faces[3 * f + j];
                    int s = faces[3 * f + (j + 1) % 3];
                    int d = faces[3 * f + (j + 2) % 3];
                    Vec3 p = (_topo.FaceCot[3 * f + j] / 3.0) * (MeshTopology.Position(_v, s) - MeshTopology.Position(_v, d));
                    Vec3 contribution = (rot[s] + rot[d] + rot[o]) * p;
                    b[s] = b[s] + contribution;
                    b[d] = b[d] - contribution;
                }
            }
            return b;
        }

        private static void ParallelFor(int count, Action<int> body)
        {
            if (count < 1000)
            {
                for (int i = 0; i < count; i++) body(i);
                return;
            }
            Parallel.ForEach(Partitioner.Create(0, count, Math.Max(256, count / (Environment.ProcessorCount * 8))), range =>
            {
                for (int i = range.Item1; i < range.Item2; i++) body(i);
            });
        }
    }
}
