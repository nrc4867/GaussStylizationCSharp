// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using GaussStylization.Internal;

namespace GaussStylization
{
    /// <summary>
    /// Gauss stylization of triangle meshes: deforms a mesh so that its face normals gather around a
    /// set of preferred directions while staying as-rigid-as-possible. Plain arrays in, plain arrays out.
    /// </summary>
    public static class GaussStylizer
    {
        /// <summary>
        /// Stylizes a triangle mesh with a single style made of discrete preferred directions
        /// (<c>Mu</c> = 1, <c>Sigma</c> = 4). For rings, other parameters or several styles use the
        /// overload that takes <see cref="GaussStyle"/> objects.
        /// </summary>
        /// <param name="vertices">Vertex positions, flat <c>[x0,y0,z0, x1,y1,z1, ...]</c> (length 3n).</param>
        /// <param name="triangles">Triangle vertex indices, flat <c>[a0,b0,c0, a1,b1,c1, ...]</c> (length 3m).
        /// The mesh must be edge-manifold (no edge shared by more than two triangles).</param>
        /// <param name="preferredDirections">The preferred surface normals, flat xyz triples (length 3k, k &gt;= 1).
        /// They are normalised internally and must be distinct. Expressed in the frame given by
        /// <paramref name="gridAxes"/> (world space when that is null).</param>
        /// <param name="vertexWeights">Optional stylization weight per vertex (length n), typically in [0,1].
        /// A face is stylized with the mean weight of its three vertices: 1 = fully, 0 = not at all
        /// (the face is then only driven by its neighbours). Null means 1 everywhere.</param>
        /// <param name="gridAxes">Optional orientation of the grid the preferred directions live in: three
        /// orthonormal axis vectors, row-major <c>[ax,ay,az, bx,by,bz, cx,cy,cz]</c> (length 9). A preferred direction
        /// <c>(d0,d1,d2)</c> then points along <c>d0*axisA + d1*axisB + d2*axisC</c> in world space. Null means the
        /// world axes (identity). Rotating the grid axes rotates the whole set of preferred directions, like
        /// the rotate keys of the original application.</param>
        /// <param name="iterations">Number of outer iterations (each is one local + global step). 0 returns a copy of the input.</param>
        /// <param name="options">Global tuning parameters; null uses the defaults.</param>
        /// <returns>New vertex positions, flat xyz (length 3n). The input arrays are not modified.</returns>
        /// <exception cref="ArgumentException">An argument has the wrong length, contains invalid values or the mesh is not edge-manifold.</exception>
        /// <exception cref="ArgumentNullException">A required array is null.</exception>
        public static double[] Stylize(
            double[] vertices,
            int[] triangles,
            double[] preferredDirections,
            double[]? vertexWeights,
            double[]? gridAxes,
            int iterations,
            GaussStylizationOptions? options = null)
        {
            if (preferredDirections == null) throw new ArgumentNullException(nameof(preferredDirections));
            return Stylize(
                vertices, triangles, new[] { new GaussStyle(preferredDirections) }, null,
                vertexWeights, gridAxes, iterations, options);
        }

        /// <summary>
        /// Stylizes a triangle mesh with one or more styles. Every triangle belongs to one style group,
        /// like the style groups painted onto the mesh in the original application.
        /// </summary>
        /// <param name="vertices">Vertex positions, flat xyz (length 3n).</param>
        /// <param name="triangles">Triangle vertex indices, flat (length 3m); the mesh must be edge-manifold.</param>
        /// <param name="styles">The styles; each triangle uses the one selected by <paramref name="faceGroups"/>.
        /// A style needs at least one direction or one ring (see <see cref="GaussStyle"/>).</param>
        /// <param name="faceGroups">Optional style index per triangle (length m, each in <c>[0, styles.Count)</c>).
        /// Null puts every triangle into style 0.</param>
        /// <param name="vertexWeights">Optional stylization weight per vertex, as in the single-style overload.</param>
        /// <param name="gridAxes">Optional grid frame, as in the single-style overload. The ring axis of every
        /// style is the grid's second axis (world Y when null).</param>
        /// <param name="iterations">Number of outer iterations.</param>
        /// <param name="options">Global tuning parameters; null uses the defaults.</param>
        /// <returns>New vertex positions, flat xyz (length 3n).</returns>
        public static double[] Stylize(
            double[] vertices,
            int[] triangles,
            IReadOnlyList<GaussStyle> styles,
            int[]? faceGroups,
            double[]? vertexWeights,
            double[]? gridAxes,
            int iterations,
            GaussStylizationOptions? options = null)
        {
            GaussStylizationCore? core = CreateCore(vertices, triangles, styles, faceGroups, vertexWeights, gridAxes, iterations, options);
            return core == null ? (double[])vertices.Clone() : core.Run(iterations);
        }

        /// <summary>
        /// Asynchronous, cancellable version of the multi-style <c>Stylize</c> overload. The argument checks run
        /// synchronously (so invalid arguments throw from this call, not from the returned task); the optimisation
        /// then runs on the thread pool. Cancellation is observed between the inner steps and inside the linear
        /// solver, so it takes effect within a fraction of one outer iteration; the returned task is then
        /// canceled (<see cref="OperationCanceledException"/>) and no partial result is returned.
        /// </summary>
        /// <param name="vertices">Vertex positions, flat xyz (length 3n).</param>
        /// <param name="triangles">Triangle vertex indices, flat (length 3m); the mesh must be edge-manifold.</param>
        /// <param name="styles">The styles; see the synchronous overload.</param>
        /// <param name="faceGroups">Optional style index per triangle; null puts every triangle into style 0.</param>
        /// <param name="vertexWeights">Optional stylization weight per vertex.</param>
        /// <param name="gridAxes">Optional grid frame (nine values).</param>
        /// <param name="iterations">Number of outer iterations.</param>
        /// <param name="options">Global tuning parameters; null uses the defaults.</param>
        /// <param name="progress">Optional; receives a report after every completed outer iteration. Reports are made
        /// from a thread-pool thread, so use <see cref="Progress{T}"/> to get them on the creating context.</param>
        /// <param name="cancellationToken">Cancels the work.</param>
        /// <returns>Task yielding new vertex positions, flat xyz (length 3n).</returns>
        public static Task<double[]> StylizeAsync(
            double[] vertices,
            int[] triangles,
            IReadOnlyList<GaussStyle> styles,
            int[]? faceGroups,
            double[]? vertexWeights,
            double[]? gridAxes,
            int iterations,
            GaussStylizationOptions? options = null,
            IProgress<GaussStylizationProgress>? progress = null,
            CancellationToken cancellationToken = default)
        {
            if (cancellationToken.IsCancellationRequested) return Task.FromCanceled<double[]>(cancellationToken);
            GaussStylizationCore? core = CreateCore(vertices, triangles, styles, faceGroups, vertexWeights, gridAxes, iterations, options);
            if (core == null) return Task.FromResult((double[])vertices.Clone());
            return Task.Run(() => core.Run(iterations, cancellationToken, progress), cancellationToken);
        }

        /// <summary>Asynchronous, cancellable version of the single-style <c>Stylize</c> overload.</summary>
        /// <param name="vertices">Vertex positions, flat xyz (length 3n).</param>
        /// <param name="triangles">Triangle vertex indices, flat (length 3m); the mesh must be edge-manifold.</param>
        /// <param name="preferredDirections">The preferred surface normals, flat xyz triples.</param>
        /// <param name="vertexWeights">Optional stylization weight per vertex.</param>
        /// <param name="gridAxes">Optional grid frame (nine values).</param>
        /// <param name="iterations">Number of outer iterations.</param>
        /// <param name="options">Global tuning parameters; null uses the defaults.</param>
        /// <param name="progress">Optional progress reports, one per completed outer iteration.</param>
        /// <param name="cancellationToken">Cancels the work.</param>
        /// <returns>Task yielding new vertex positions, flat xyz (length 3n).</returns>
        public static Task<double[]> StylizeAsync(
            double[] vertices,
            int[] triangles,
            double[] preferredDirections,
            double[]? vertexWeights,
            double[]? gridAxes,
            int iterations,
            GaussStylizationOptions? options = null,
            IProgress<GaussStylizationProgress>? progress = null,
            CancellationToken cancellationToken = default)
        {
            if (preferredDirections == null) throw new ArgumentNullException(nameof(preferredDirections));
            return StylizeAsync(
                vertices, triangles, new[] { new GaussStyle(preferredDirections) }, null,
                vertexWeights, gridAxes, iterations, options, progress, cancellationToken);
        }

        /// <summary>Validates everything; returns null when the result is simply a copy of the input.</summary>
        private static GaussStylizationCore? CreateCore(
            double[] vertices,
            int[] triangles,
            IReadOnlyList<GaussStyle> styles,
            int[]? faceGroups,
            double[]? vertexWeights,
            double[]? gridAxes,
            int iterations,
            GaussStylizationOptions? options)
        {
            if (vertices == null) throw new ArgumentNullException(nameof(vertices));
            if (triangles == null) throw new ArgumentNullException(nameof(triangles));
            if (styles == null) throw new ArgumentNullException(nameof(styles));
            options ??= new GaussStylizationOptions();

            ValidateOptions(options);
            if (iterations < 0) throw new ArgumentOutOfRangeException(nameof(iterations), "iterations must be >= 0.");
            if (vertices.Length % 3 != 0) throw new ArgumentException("vertices length must be a multiple of 3.", nameof(vertices));
            if (triangles.Length % 3 != 0) throw new ArgumentException("triangles length must be a multiple of 3.", nameof(triangles));
            RequireFinite(vertices, nameof(vertices));

            int nv = vertices.Length / 3;
            int nf = triangles.Length / 3;
            foreach (int idx in triangles)
                if (idx < 0 || idx >= nv) throw new ArgumentException("triangles contains a vertex index out of range.", nameof(triangles));
            for (int f = 0; f < nf; f++)
            {
                int a = triangles[3 * f], b = triangles[3 * f + 1], c = triangles[3 * f + 2];
                if (a == b || b == c || a == c)
                    throw new ArgumentException($"Triangle {f} references the same vertex more than once.", nameof(triangles));
            }

            if (styles.Count == 0) throw new ArgumentException("At least one style is required.", nameof(styles));

            if (faceGroups != null)
            {
                if (faceGroups.Length != nf) throw new ArgumentException("faceGroups must have one entry per triangle.", nameof(faceGroups));
                foreach (int g in faceGroups)
                    if (g < 0 || g >= styles.Count) throw new ArgumentException("faceGroups contains a style index out of range.", nameof(faceGroups));
            }

            if (vertexWeights != null)
            {
                if (vertexWeights.Length != nv) throw new ArgumentException("vertexWeights must have one entry per vertex.", nameof(vertexWeights));
                foreach (double w in vertexWeights)
                    if (double.IsNaN(w) || double.IsInfinity(w) || w < 0)
                        throw new ArgumentException("vertexWeights must be finite and >= 0.", nameof(vertexWeights));
            }

            Mat3 grid = Mat3.Identity;
            if (gridAxes != null)
            {
                if (gridAxes.Length != 9) throw new ArgumentException("gridAxes must contain 9 values (three axis vectors).", nameof(gridAxes));
                RequireFinite(gridAxes, nameof(gridAxes));
                var a0 = new Vec3(gridAxes[0], gridAxes[1], gridAxes[2]);
                var a1 = new Vec3(gridAxes[3], gridAxes[4], gridAxes[5]);
                var a2 = new Vec3(gridAxes[6], gridAxes[7], gridAxes[8]);
                const double tol = 1e-4;
                if (Math.Abs(a0.Norm - 1) > tol || Math.Abs(a1.Norm - 1) > tol || Math.Abs(a2.Norm - 1) > tol ||
                    Math.Abs(a0.Dot(a1)) > tol || Math.Abs(a0.Dot(a2)) > tol || Math.Abs(a1.Dot(a2)) > tol)
                    throw new ArgumentException("gridAxes must be three orthonormal vectors.", nameof(gridAxes));
                grid = Mat3.FromColumns(a0, a1, a2); // world = A * local, columns are the axes
            }
            Vec3 ringAxis = (grid * new Vec3(0, 1, 0)).Normalized();

            if (options.FixedVertices != null)
                foreach (int v in options.FixedVertices)
                    if (v < 0 || v >= nv) throw new ArgumentException("FixedVertices contains an index out of range.", nameof(options));

            var functions = new GFunction[styles.Count];
            var mu = new double[styles.Count];
            for (int s = 0; s < styles.Count; s++)
            {
                GaussStyle style = styles[s] ?? throw new ArgumentException($"Style {s} is null.", nameof(styles));
                functions[s] = BuildFunction(style, s, grid, ringAxis);
                mu[s] = style.Mu;
            }

            if (iterations == 0 || nf == 0) return null;

            var topoProbe = new MeshTopology(vertices, triangles); // also validates manifoldness
            bool[] fixedVertices = topoProbe.ComputeFixedVertices(options.FixedVertices);

            int[] groups = faceGroups ?? new int[nf];
            return new GaussStylizationCore(vertices, triangles, functions, mu, groups, vertexWeights, options, fixedVertices);
        }

        /// <summary>Single-precision convenience overload; computation is done in double precision.</summary>
        public static float[] Stylize(
            float[] vertices,
            int[] triangles,
            float[] preferredDirections,
            float[]? vertexWeights,
            float[]? gridAxes,
            int iterations,
            GaussStylizationOptions? options = null)
        {
            if (vertices == null) throw new ArgumentNullException(nameof(vertices));
            if (preferredDirections == null) throw new ArgumentNullException(nameof(preferredDirections));
            double[] result = Stylize(
                ToDouble(vertices)!, triangles, ToDouble(preferredDirections)!,
                ToDouble(vertexWeights), ToDouble(gridAxes), iterations, options);
            var output = new float[result.Length];
            for (int i = 0; i < output.Length; i++) output[i] = (float)result[i];
            return output;
        }

        /// <summary>The six axis-aligned cube normals (+-X, +-Y, +-Z), the default style of the original application.</summary>
        public static double[] CubeDirections() => new double[]
        {
            0, 0, 1,  0, 0, -1,
            0, 1, 0,  0, -1, 0,
            1, 0, 0,  -1, 0, 0,
        };

        private static GFunction BuildFunction(GaussStyle style, int index, Mat3 grid, Vec3 ringAxis)
        {
            string name = $"styles[{index}]";
            double[] dirs = style.Directions ?? Array.Empty<double>();
            double[] rings = style.RingHeights ?? Array.Empty<double>();

            if (dirs.Length % 3 != 0) throw new ArgumentException($"{name}.Directions must contain xyz triples.", nameof(style));
            RequireFinite(dirs, name + ".Directions");
            RequireFinite(rings, name + ".RingHeights");
            if (dirs.Length == 0 && rings.Length == 0)
                throw new ArgumentException($"{name} needs at least one direction or one ring.", nameof(style));
            foreach (double d in rings)
                if (d < -1 || d > 1) throw new ArgumentException($"{name}.RingHeights must lie in [-1, 1].", nameof(style));
            if (!(style.Mu >= 0) || double.IsInfinity(style.Mu)) throw new ArgumentException($"{name}.Mu must be finite and >= 0.", nameof(style));
            if (!(style.Sigma > 0) || double.IsInfinity(style.Sigma)) throw new ArgumentException($"{name}.Sigma must be finite and > 0.", nameof(style));
            if (!(style.RingContribution >= 0) || double.IsInfinity(style.RingContribution))
                throw new ArgumentException($"{name}.RingContribution must be finite and >= 0.", nameof(style));

            var normals = new Vec3[dirs.Length / 3];
            for (int k = 0; k < normals.Length; k++)
            {
                var local = new Vec3(dirs[3 * k], dirs[3 * k + 1], dirs[3 * k + 2]);
                if (local.Norm == 0) throw new ArgumentException($"{name}: direction {k} has zero length.", nameof(style));
                normals[k] = (grid * local).Normalized();
            }
            return new GFunction(normals, style.Sigma, (double[])rings.Clone(), ringAxis, style.RingContribution);
        }

        private static double[]? ToDouble(float[]? a)
        {
            if (a == null) return null;
            var r = new double[a.Length];
            for (int i = 0; i < r.Length; i++) r[i] = a[i];
            return r;
        }

        private static void RequireFinite(double[] a, string name)
        {
            foreach (double x in a)
                if (double.IsNaN(x) || double.IsInfinity(x))
                    throw new ArgumentException("Values must be finite.", name);
        }

        private static void ValidateOptions(GaussStylizationOptions o)
        {
            if (!(o.Lambda >= 0) || double.IsInfinity(o.Lambda)) throw new ArgumentException("Lambda must be finite and >= 0.", nameof(o));
            if (o.InnerIterations < 1) throw new ArgumentException("InnerIterations must be >= 1.", nameof(o));
            if (!(o.SolverTolerance > 0)) throw new ArgumentException("SolverTolerance must be > 0.", nameof(o));
            if (o.SolverMaxIterations < 1) throw new ArgumentException("SolverMaxIterations must be >= 1.", nameof(o));
        }
    }
}
