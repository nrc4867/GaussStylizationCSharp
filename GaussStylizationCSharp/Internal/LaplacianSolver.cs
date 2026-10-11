// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

using System;
using System.Threading;

namespace GaussStylization.Internal
{
    /// <summary>
    /// Solves the global-step system of gaussStylization, <c>L U = -B</c> with some vertices fixed
    /// (igl::min_quad_with_fixed with A = cotangent Laplacian, no equality constraints), using
    /// Jacobi-preconditioned conjugate gradients on the reduced, positive definite system -L_ff.
    /// This replaces Eigen's sparse direct factorisation with plain managed code.
    /// </summary>
    internal sealed class LaplacianSolver
    {
        private readonly int _n;
        private readonly int[] _unknownIndex;     // vertex -> reduced index, -1 if fixed
        private readonly int[] _unknownVertex;    // reduced index -> vertex
        private readonly int[] _rowStart;
        private readonly int[] _col;
        private readonly double[] _val;
        private readonly double[] _diag;
        private readonly int[] _fixedRowStart;    // for each unknown: fixed neighbours (for the RHS)
        private readonly int[] _fixedVertex;
        private readonly double[] _fixedWeight;

        private readonly double _tolerance;
        private readonly int _maxIterations;

        public LaplacianSolver(MeshTopology topo, bool[] isFixed, double tolerance, int maxIterations)
        {
            _n = topo.VertexCount;
            _tolerance = tolerance;
            _maxIterations = maxIterations;

            _unknownIndex = new int[_n];
            int nu = 0;
            for (int v = 0; v < _n; v++) _unknownIndex[v] = isFixed[v] ? -1 : nu++;
            _unknownVertex = new int[nu];
            for (int v = 0; v < _n; v++) if (_unknownIndex[v] >= 0) _unknownVertex[_unknownIndex[v]] = v;

            var offCount = new int[nu + 1];
            var fixedCount = new int[nu + 1];
            _diag = new double[nu];
            for (int e = 0; e < topo.EdgeCount; e++)
            {
                int a = topo.EdgeVertices[2 * e], b = topo.EdgeVertices[2 * e + 1];
                int ia = _unknownIndex[a], ib = _unknownIndex[b];
                double w = topo.EdgeWeight[e];
                if (ia >= 0) { _diag[ia] += w; if (ib >= 0) offCount[ia + 1]++; else fixedCount[ia + 1]++; }
                if (ib >= 0) { _diag[ib] += w; if (ia >= 0) offCount[ib + 1]++; else fixedCount[ib + 1]++; }
            }
            for (int i = 0; i < nu; i++) { offCount[i + 1] += offCount[i]; fixedCount[i + 1] += fixedCount[i]; }
            _rowStart = offCount;
            _fixedRowStart = fixedCount;
            _col = new int[_rowStart[nu]];
            _val = new double[_rowStart[nu]];
            _fixedVertex = new int[_fixedRowStart[nu]];
            _fixedWeight = new double[_fixedRowStart[nu]];

            var oc = new int[nu];
            var fc = new int[nu];
            for (int e = 0; e < topo.EdgeCount; e++)
            {
                int a = topo.EdgeVertices[2 * e], b = topo.EdgeVertices[2 * e + 1];
                int ia = _unknownIndex[a], ib = _unknownIndex[b];
                double w = topo.EdgeWeight[e];
                if (ia >= 0)
                {
                    if (ib >= 0) { int k = _rowStart[ia] + oc[ia]++; _col[k] = ib; _val[k] = -w; }
                    else { int k = _fixedRowStart[ia] + fc[ia]++; _fixedVertex[k] = b; _fixedWeight[k] = w; }
                }
                if (ib >= 0)
                {
                    if (ia >= 0) { int k = _rowStart[ib] + oc[ib]++; _col[k] = ia; _val[k] = -w; }
                    else { int k = _fixedRowStart[ib] + fc[ib]++; _fixedVertex[k] = a; _fixedWeight[k] = w; }
                }
            }
        }

        /// <summary>
        /// Solves for one coordinate. <paramref name="b"/> is the per-vertex right-hand side B,
        /// <paramref name="fixedValues"/> the prescribed values (read at fixed vertices), and
        /// <paramref name="x"/> holds the initial guess on entry and the full solution on exit.
        /// </summary>
        public void Solve(double[] b, double[] fixedValues, double[] x, CancellationToken cancellationToken = default)
        {
            int nu = _unknownVertex.Length;
            var rhs = new double[nu];
            var sol = new double[nu];
            var r = new double[nu];
            var z = new double[nu];
            var p = new double[nu];
            var q = new double[nu];
            var invDiag = new double[nu];

            for (int i = 0; i < nu; i++)
            {
                double s = b[_unknownVertex[i]];
                for (int k = _fixedRowStart[i]; k < _fixedRowStart[i + 1]; k++)
                    s += _fixedWeight[k] * fixedValues[_fixedVertex[k]];
                rhs[i] = s;
                sol[i] = x[_unknownVertex[i]];
                invDiag[i] = _diag[i] > 0 ? 1.0 / _diag[i] : 1.0;
            }

            double bnorm = Math.Sqrt(Dot(rhs, rhs));
            if (bnorm > 0)
            {
                Multiply(sol, q);
                for (int i = 0; i < nu; i++) r[i] = rhs[i] - q[i];
                double target = _tolerance * bnorm;
                if (Math.Sqrt(Dot(r, r)) > target)
                {
                    for (int i = 0; i < nu; i++) { z[i] = invDiag[i] * r[i]; p[i] = z[i]; }
                    double rz = Dot(r, z);
                    for (int it = 0; it < _maxIterations; it++)
                    {
                        cancellationToken.ThrowIfCancellationRequested();
                        Multiply(p, q);
                        double pq = Dot(p, q);
                        if (!(pq > 0)) break;
                        double alpha = rz / pq;
                        for (int i = 0; i < nu; i++) { sol[i] += alpha * p[i]; r[i] -= alpha * q[i]; }
                        if (Math.Sqrt(Dot(r, r)) <= target) break;
                        for (int i = 0; i < nu; i++) z[i] = invDiag[i] * r[i];
                        double rzNew = Dot(r, z);
                        double beta = rzNew / rz;
                        rz = rzNew;
                        for (int i = 0; i < nu; i++) p[i] = z[i] + beta * p[i];
                    }
                }
            }

            for (int v = 0; v < _n; v++)
                if (_unknownIndex[v] < 0) x[v] = fixedValues[v];
            for (int i = 0; i < nu; i++) x[_unknownVertex[i]] = sol[i];
        }

        private void Multiply(double[] v, double[] result)
        {
            for (int i = 0; i < result.Length; i++)
            {
                double s = _diag[i] * v[i];
                for (int k = _rowStart[i]; k < _rowStart[i + 1]; k++) s += _val[k] * v[_col[k]];
                result[i] = s;
            }
        }

        private static double Dot(double[] a, double[] b)
        {
            double s = 0;
            for (int i = 0; i < a.Length; i++) s += a[i] * b[i];
            return s;
        }
    }
}
