// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

namespace GaussStylizationCSharp.Tests;

internal static class TestMeshes
{
    /// <summary>Unit icosphere with the given number of subdivisions (20 * 4^s faces).</summary>
    public static (double[] V, int[] F) IcoSphere(int subdivisions)
    {
        double t = (1.0 + Math.Sqrt(5.0)) / 2.0;
        var verts = new List<double[]>
        {
            new[] { -1, t, 0.0 }, new[] { 1, t, 0.0 }, new[] { -1, -t, 0.0 }, new[] { 1, -t, 0.0 },
            new[] { 0.0, -1, t }, new[] { 0.0, 1, t }, new[] { 0.0, -1, -t }, new[] { 0.0, 1, -t },
            new[] { t, 0.0, -1 }, new[] { t, 0.0, 1 }, new[] { -t, 0.0, -1 }, new[] { -t, 0.0, 1 },
        };
        for (int i = 0; i < verts.Count; i++) verts[i] = Normalize(verts[i]);
        var faces = new List<int[]>
        {
            new[] { 0, 11, 5 }, new[] { 0, 5, 1 }, new[] { 0, 1, 7 }, new[] { 0, 7, 10 }, new[] { 0, 10, 11 },
            new[] { 1, 5, 9 }, new[] { 5, 11, 4 }, new[] { 11, 10, 2 }, new[] { 10, 7, 6 }, new[] { 7, 1, 8 },
            new[] { 3, 9, 4 }, new[] { 3, 4, 2 }, new[] { 3, 2, 6 }, new[] { 3, 6, 8 }, new[] { 3, 8, 9 },
            new[] { 4, 9, 5 }, new[] { 2, 4, 11 }, new[] { 6, 2, 10 }, new[] { 8, 6, 7 }, new[] { 9, 8, 1 },
        };
        for (int s = 0; s < subdivisions; s++)
        {
            var cache = new Dictionary<(int, int), int>();
            int Mid(int a, int b)
            {
                var key = a < b ? (a, b) : (b, a);
                if (cache.TryGetValue(key, out int m)) return m;
                var p = Normalize(new[]
                {
                    (verts[a][0] + verts[b][0]) / 2, (verts[a][1] + verts[b][1]) / 2, (verts[a][2] + verts[b][2]) / 2,
                });
                verts.Add(p);
                return cache[key] = verts.Count - 1;
            }
            var next = new List<int[]>();
            foreach (var f in faces)
            {
                int ab = Mid(f[0], f[1]), bc = Mid(f[1], f[2]), ca = Mid(f[2], f[0]);
                next.Add(new[] { f[0], ab, ca });
                next.Add(new[] { f[1], bc, ab });
                next.Add(new[] { f[2], ca, bc });
                next.Add(new[] { ab, bc, ca });
            }
            faces = next;
        }
        return (verts.SelectMany(v => v).ToArray(), faces.SelectMany(f => f).ToArray());
    }

    /// <summary>Flat n x n grid of quads (two triangles each) in the z = 0 plane, spacing 1.</summary>
    public static (double[] V, int[] F) Grid(int n)
    {
        var v = new double[3 * n * n];
        for (int y = 0; y < n; y++)
            for (int x = 0; x < n; x++)
            {
                int i = y * n + x;
                v[3 * i] = x; v[3 * i + 1] = y; v[3 * i + 2] = 0;
            }
        var f = new List<int>();
        for (int y = 0; y < n - 1; y++)
            for (int x = 0; x < n - 1; x++)
            {
                int a = y * n + x, b = a + 1, c = a + n, d = c + 1;
                f.AddRange(new[] { a, b, d, a, d, c });
            }
        return (v, f.ToArray());
    }

    public static double[] Normalize(double[] p)
    {
        double l = Math.Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        return new[] { p[0] / l, p[1] / l, p[2] / l };
    }

    /// <summary>Mean over faces of max_k (n_f . d_k): how well face normals gather around the directions.</summary>
    public static double Alignment(double[] v, int[] f, double[] dirs, Func<int, bool>? faceFilter = null)
    {
        double sum = 0;
        int nf = f.Length / 3, counted = 0;
        for (int i = 0; i < nf; i++)
        {
            if (faceFilter != null && !faceFilter(i)) continue;
            counted++;
            double[] p0 = P(v, f[3 * i]), p1 = P(v, f[3 * i + 1]), p2 = P(v, f[3 * i + 2]);
            double[] e1 = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
            double[] e2 = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
            var n = Normalize(new[]
            {
                e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0],
            });
            double best = -1;
            for (int k = 0; k < dirs.Length / 3; k++)
            {
                var d = Normalize(new[] { dirs[3 * k], dirs[3 * k + 1], dirs[3 * k + 2] });
                best = Math.Max(best, n[0] * d[0] + n[1] * d[1] + n[2] * d[2]);
            }
            sum += best;
        }
        return sum / counted;
    }

    public static double MaxAbsDiff(double[] a, double[] b)
    {
        double m = 0;
        for (int i = 0; i < a.Length; i++) m = Math.Max(m, Math.Abs(a[i] - b[i]));
        return m;
    }

    private static double[] P(double[] v, int i) => new[] { v[3 * i], v[3 * i + 1], v[3 * i + 2] };
}
