// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.
// Port of utils_gauss/g_function.cpp and utils_gauss/normalize_g.cpp.

using System;

namespace GaussStylization.Internal
{
    /// <summary>
    /// The "g function" of Gauss stylization. The discrete part is
    /// <c>sum_k w_k exp(sigma n_k . x)</c> with weights chosen so that g(n_k) = 1 at every preferred
    /// normal n_k. The optional semi-discrete part adds rings around an axis a:
    /// <c>g(x) = c * discrete(x) + sum_r v_r exp(sigma (1 - (a.x - d_r)^2))</c>, where d_r are the ring
    /// heights, c the ring contribution, and the v_r are chosen so the ring terms equal 1 at each ring.
    /// As in the original, <c>c</c> only applies when at least one ring exists.
    /// </summary>
    internal sealed class GFunction
    {
        private readonly Vec3[] _normals;
        private readonly double[] _weights;
        private readonly double _sigma;

        private readonly double[] _ringHeights;
        private readonly double[] _ringWeights;
        private readonly Vec3 _axis;
        private readonly double _ringContribution;

        public GFunction(Vec3[] unitNormals, double sigma)
            : this(unitNormals, sigma, Array.Empty<double>(), new Vec3(0, 1, 0), 1.0)
        {
        }

        public GFunction(Vec3[] unitNormals, double sigma, double[] ringHeights, Vec3 unitAxis, double ringContribution)
        {
            _normals = unitNormals;
            _sigma = sigma;
            _ringHeights = ringHeights;
            _axis = unitAxis;
            _ringContribution = ringContribution;

            // normalize_g, discrete part: solve A w = 1 with A_ij = exp(sigma n_i . n_j)
            int k = unitNormals.Length;
            var a = new double[k, k + 1];
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < k; j++) a[i, j] = Math.Exp(sigma * unitNormals[i].Dot(unitNormals[j]));
                a[i, k] = 1.0;
            }
            _weights = SolveDense(a, k);

            // normalize_g, ring part: solve B v = 1 with B_ij = exp(sigma (1 - (d_i - d_j)^2))
            int r = ringHeights.Length;
            var b = new double[r, r + 1];
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < r; j++)
                {
                    double dd = ringHeights[i] - ringHeights[j];
                    b[i, j] = Math.Exp(sigma * (1.0 - dd * dd));
                }
                b[i, r] = 1.0;
            }
            _ringWeights = SolveDense(b, r);
        }

        public double Value(Vec3 x)
        {
            double v = 0;
            for (int k = 0; k < _normals.Length; k++) v += _weights[k] * Math.Exp(_sigma * _normals[k].Dot(x));
            if (_ringHeights.Length > 0)
            {
                v *= _ringContribution;
                double aw = _axis.Dot(x);
                for (int k = 0; k < _ringHeights.Length; k++)
                {
                    double t = aw - _ringHeights[k];
                    v += _ringWeights[k] * Math.Exp(_sigma * (1.0 - t * t));
                }
            }
            return v;
        }

        /// <summary>Gradient and Hessian of g at x (equations 8 and 9 of the paper).</summary>
        public void GradientAndHessian(Vec3 x, out Vec3 grad, out Mat3 hess)
        {
            grad = new Vec3(0, 0, 0);
            hess = default;
            for (int k = 0; k < _normals.Length; k++)
            {
                Vec3 n = _normals[k];
                double e = _weights[k] * Math.Exp(_sigma * n.Dot(x));
                grad = grad + (_sigma * e) * n;
                hess = hess + (_sigma * _sigma * e) * Mat3.Outer(n, n);
            }

            if (_ringHeights.Length > 0)
            {
                grad = _ringContribution * grad;
                hess = _ringContribution * hess;
                double aw = _axis.Dot(x);
                Mat3 aa = Mat3.Outer(_axis, _axis);
                for (int k = 0; k < _ringHeights.Length; k++)
                {
                    double t = aw - _ringHeights[k];
                    double e = _ringWeights[k] * Math.Exp(_sigma * (1.0 - t * t));
                    grad = grad + (-2.0 * _sigma * t * e) * _axis;
                    hess = hess + (2.0 * _sigma * (2.0 * _sigma * t * t - 1.0) * e) * aa;
                }
            }
        }

        // Gaussian elimination with partial pivoting on the augmented matrix [A | b].
        private static double[] SolveDense(double[,] m, int n)
        {
            double scale = 0;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) scale = Math.Max(scale, Math.Abs(m[i, j]));
            for (int c = 0; c < n; c++)
            {
                int piv = c;
                for (int r = c + 1; r < n; r++) if (Math.Abs(m[r, c]) > Math.Abs(m[piv, c])) piv = r;
                if (!(Math.Abs(m[piv, c]) > 1e-12 * scale))
                    throw new ArgumentException(
                        "The preferred directions or ring heights are (numerically) linearly dependent for this sigma; " +
                        "remove duplicate or near-duplicate entries.");
                if (piv != c)
                    for (int j = c; j <= n; j++) { double t = m[c, j]; m[c, j] = m[piv, j]; m[piv, j] = t; }
                for (int r = c + 1; r < n; r++)
                {
                    double f = m[r, c] / m[c, c];
                    if (f == 0) continue;
                    for (int j = c; j <= n; j++) m[r, j] -= f * m[c, j];
                }
            }
            var x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                double s = m[i, n];
                for (int j = i + 1; j < n; j++) s -= m[i, j] * x[j];
                x[i] = s / m[i, i];
                if (double.IsNaN(x[i]) || double.IsInfinity(x[i]))
                    throw new ArgumentException("Could not normalise the g function; check sigma, the directions and the ring heights.");
            }
            return x;
        }
    }
}
