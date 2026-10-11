// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

using System;

namespace GaussStylization.Internal
{
    /// <summary>Row-major 3x3 matrix.</summary>
    internal struct Mat3
    {
        public double M00, M01, M02, M10, M11, M12, M20, M21, M22;

        public static Mat3 Identity => new Mat3 { M00 = 1, M11 = 1, M22 = 1 };

        public static Mat3 FromColumns(Vec3 c0, Vec3 c1, Vec3 c2) => new Mat3
        {
            M00 = c0.X, M01 = c1.X, M02 = c2.X,
            M10 = c0.Y, M11 = c1.Y, M12 = c2.Y,
            M20 = c0.Z, M21 = c1.Z, M22 = c2.Z,
        };

        public static Mat3 Outer(Vec3 a, Vec3 b) => new Mat3
        {
            M00 = a.X * b.X, M01 = a.X * b.Y, M02 = a.X * b.Z,
            M10 = a.Y * b.X, M11 = a.Y * b.Y, M12 = a.Y * b.Z,
            M20 = a.Z * b.X, M21 = a.Z * b.Y, M22 = a.Z * b.Z,
        };

        public static Mat3 operator +(Mat3 a, Mat3 b) => new Mat3
        {
            M00 = a.M00 + b.M00, M01 = a.M01 + b.M01, M02 = a.M02 + b.M02,
            M10 = a.M10 + b.M10, M11 = a.M11 + b.M11, M12 = a.M12 + b.M12,
            M20 = a.M20 + b.M20, M21 = a.M21 + b.M21, M22 = a.M22 + b.M22,
        };

        public static Mat3 operator -(Mat3 a, Mat3 b) => new Mat3
        {
            M00 = a.M00 - b.M00, M01 = a.M01 - b.M01, M02 = a.M02 - b.M02,
            M10 = a.M10 - b.M10, M11 = a.M11 - b.M11, M12 = a.M12 - b.M12,
            M20 = a.M20 - b.M20, M21 = a.M21 - b.M21, M22 = a.M22 - b.M22,
        };

        public static Mat3 operator *(double s, Mat3 a) => new Mat3
        {
            M00 = s * a.M00, M01 = s * a.M01, M02 = s * a.M02,
            M10 = s * a.M10, M11 = s * a.M11, M12 = s * a.M12,
            M20 = s * a.M20, M21 = s * a.M21, M22 = s * a.M22,
        };

        public static Vec3 operator *(Mat3 a, Vec3 v) => new Vec3(
            a.M00 * v.X + a.M01 * v.Y + a.M02 * v.Z,
            a.M10 * v.X + a.M11 * v.Y + a.M12 * v.Z,
            a.M20 * v.X + a.M21 * v.Y + a.M22 * v.Z);

        public static Mat3 operator *(Mat3 a, Mat3 b) => new Mat3
        {
            M00 = a.M00 * b.M00 + a.M01 * b.M10 + a.M02 * b.M20,
            M01 = a.M00 * b.M01 + a.M01 * b.M11 + a.M02 * b.M21,
            M02 = a.M00 * b.M02 + a.M01 * b.M12 + a.M02 * b.M22,
            M10 = a.M10 * b.M00 + a.M11 * b.M10 + a.M12 * b.M20,
            M11 = a.M10 * b.M01 + a.M11 * b.M11 + a.M12 * b.M21,
            M12 = a.M10 * b.M02 + a.M11 * b.M12 + a.M12 * b.M22,
            M20 = a.M20 * b.M00 + a.M21 * b.M10 + a.M22 * b.M20,
            M21 = a.M20 * b.M01 + a.M21 * b.M11 + a.M22 * b.M21,
            M22 = a.M20 * b.M02 + a.M21 * b.M12 + a.M22 * b.M22,
        };

        public Mat3 Transposed() => new Mat3
        {
            M00 = M00, M01 = M10, M02 = M20,
            M10 = M01, M11 = M11, M12 = M21,
            M20 = M02, M21 = M12, M22 = M22,
        };

        public double Determinant =>
            M00 * (M11 * M22 - M12 * M21) - M01 * (M10 * M22 - M12 * M20) + M02 * (M10 * M21 - M11 * M20);

        public Vec3 Column(int j) =>
            j == 0 ? new Vec3(M00, M10, M20) : (j == 1 ? new Vec3(M01, M11, M21) : new Vec3(M02, M12, M22));

        /// <summary>
        /// Returns the proper rotation R (det +1) maximising trace(R * S). For S = U Sigma V^T this is
        /// R = V diag(1,1,det(V U^T)) U^T, i.e. the transpose of what <c>igl::polar_svd3x3</c> returns
        /// (which gaussStylization transposes right after the call). A one-sided Jacobi SVD is used,
        /// which stays accurate for rank-deficient input.
        /// </summary>
        public static Mat3 FitRotation(Mat3 s)
        {
            Vec3 a0 = s.Column(0), a1 = s.Column(1), a2 = s.Column(2);
            Vec3 v0 = new Vec3(1, 0, 0), v1 = new Vec3(0, 1, 0), v2 = new Vec3(0, 0, 1);

            for (int sweep = 0; sweep < 40; sweep++)
            {
                bool rotated = false;
                JacobiRotate(ref a0, ref a1, ref v0, ref v1, ref rotated);
                JacobiRotate(ref a0, ref a2, ref v0, ref v2, ref rotated);
                JacobiRotate(ref a1, ref a2, ref v1, ref v2, ref rotated);
                if (!rotated) break;
            }

            // Sort singular values (column norms) in descending order.
            double s0 = a0.Norm, s1 = a1.Norm, s2 = a2.Norm;
            if (s0 < s1) { Swap(ref s0, ref s1); Swap(ref a0, ref a1); Swap(ref v0, ref v1); }
            if (s0 < s2) { Swap(ref s0, ref s2); Swap(ref a0, ref a2); Swap(ref v0, ref v2); }
            if (s1 < s2) { Swap(ref s1, ref s2); Swap(ref a1, ref a2); Swap(ref v1, ref v2); }

            if (!(s0 > 1e-300)) return Identity;

            Vec3 u0 = a0 * (1.0 / s0);
            Vec3 u1;
            if (s1 > 1e-12 * s0)
            {
                u1 = (a1 - u0.Dot(a1) * u0).Normalized();
            }
            else
            {
                // Rank one: any unit vector orthogonal to u0 will do.
                Vec3 t = Math.Abs(u0.X) < 0.9 ? new Vec3(1, 0, 0) : new Vec3(0, 1, 0);
                u1 = (t - u0.Dot(t) * u0).Normalized();
            }
            Vec3 u2 = u0.Cross(u1); // det(U) = +1

            double d = FromColumns(v0, v1, v2).Determinant >= 0 ? 1.0 : -1.0;
            Mat3 vd = FromColumns(v0, v1, d * v2);
            return vd * FromColumns(u0, u1, u2).Transposed();
        }

        /// <summary>
        /// Solves A x = b for symmetric A with an LDL^T factorisation using diagonal pivoting.
        /// A numerically zero pivot yields a zero component (as Eigen's LDLT does), so singular
        /// systems return a finite least-effort solution instead of NaN.
        /// </summary>
        public static Vec3 SolveSymmetric(Mat3 m, Vec3 rhs)
        {
            Span<double> a = stackalloc double[9];
            a[0] = m.M00; a[1] = m.M01; a[2] = m.M02;
            a[3] = m.M10; a[4] = m.M11; a[5] = m.M12;
            a[6] = m.M20; a[7] = m.M21; a[8] = m.M22;
            Span<double> b = stackalloc double[3];
            b[0] = rhs.X; b[1] = rhs.Y; b[2] = rhs.Z;
            Span<int> p = stackalloc int[3];
            p[0] = 0; p[1] = 1; p[2] = 2;
            Span<double> dg = stackalloc double[3];
            Span<double> l = stackalloc double[9];
            l.Clear();

            double scale = Math.Max(Math.Abs(a[0]), Math.Max(Math.Abs(a[4]), Math.Abs(a[8])));
            double tol = scale * 1e-13;

            for (int k = 0; k < 3; k++)
            {
                int piv = k;
                for (int i = k + 1; i < 3; i++)
                    if (Math.Abs(a[i * 3 + i]) > Math.Abs(a[piv * 3 + piv])) piv = i;
                if (piv != k)
                {
                    for (int j = 0; j < 3; j++) Swap(ref a[k * 3 + j], ref a[piv * 3 + j]);
                    for (int i = 0; i < 3; i++) Swap(ref a[i * 3 + k], ref a[i * 3 + piv]);
                    for (int j = 0; j < k; j++) Swap(ref l[k * 3 + j], ref l[piv * 3 + j]);
                    Swap(ref b[k], ref b[piv]);
                    Swap(ref p[k], ref p[piv]);
                }
                double d = a[k * 3 + k];
                if (Math.Abs(d) <= tol) { dg[k] = 0; continue; }
                dg[k] = d;
                for (int i = k + 1; i < 3; i++) l[i * 3 + k] = a[i * 3 + k] / d;
                for (int i = k + 1; i < 3; i++)
                    for (int j = k + 1; j < 3; j++)
                        a[i * 3 + j] -= l[i * 3 + k] * a[k * 3 + j];
            }

            Span<double> y = stackalloc double[3];
            for (int i = 0; i < 3; i++)
            {
                double sum = b[i];
                for (int k = 0; k < i; k++) sum -= l[i * 3 + k] * y[k];
                y[i] = sum;
            }
            for (int i = 0; i < 3; i++) y[i] = dg[i] != 0 ? y[i] / dg[i] : 0.0;
            Span<double> x = stackalloc double[3];
            for (int i = 2; i >= 0; i--)
            {
                double sum = y[i];
                for (int k = i + 1; k < 3; k++) sum -= l[k * 3 + i] * x[k];
                x[i] = sum;
            }
            Span<double> res = stackalloc double[3];
            for (int i = 0; i < 3; i++) res[p[i]] = x[i];
            return new Vec3(res[0], res[1], res[2]);
        }

        private static void Swap<T>(ref T a, ref T b) { T t = a; a = b; b = t; }

        private static void JacobiRotate(ref Vec3 ap, ref Vec3 aq, ref Vec3 vp, ref Vec3 vq, ref bool rotated)
        {
            double alpha = ap.SquaredNorm, beta = aq.SquaredNorm, gamma = ap.Dot(aq);
            if (gamma == 0 || Math.Abs(gamma) <= 1e-15 * Math.Sqrt(alpha * beta)) return;
            double zeta = (beta - alpha) / (2.0 * gamma);
            double t = (zeta >= 0 ? 1.0 : -1.0) / (Math.Abs(zeta) + Math.Sqrt(1.0 + zeta * zeta));
            double c = 1.0 / Math.Sqrt(1.0 + t * t);
            double sn = c * t;
            Vec3 nap = c * ap - sn * aq, naq = sn * ap + c * aq;
            Vec3 nvp = c * vp - sn * vq, nvq = sn * vp + c * vq;
            ap = nap; aq = naq; vp = nvp; vq = nvq;
            rotated = true;
        }
    }
}
