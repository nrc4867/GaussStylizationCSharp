// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

namespace GaussStylization
{
    /// <summary>
    /// Global tuning parameters of the stylization. The defaults are those used by the original
    /// gaussStylization application (lambda = 4, one inner iteration). Per-style parameters
    /// (<c>Mu</c>, <c>Sigma</c>, ring weight) live on <see cref="GaussStyle"/>.
    /// </summary>
    public sealed class GaussStylizationOptions
    {
        /// <summary>
        /// Coupling between the rotated-edge (ARAP) term and the stylization target. Larger values
        /// keep the surface closer to the target normals / edges. Must be &gt;= 0. Default 4.
        /// </summary>
        public double Lambda { get; set; } = 4.0;

        /// <summary>
        /// Number of face/edge (ADMM) updates performed inside each outer iteration. The original
        /// application defaults to 1; its Python example uses 6. Must be &gt;= 1.
        /// </summary>
        public int InnerIterations { get; set; } = 1;

        /// <summary>
        /// Indices of vertices that stay at their input position. If null or empty, one vertex per
        /// connected component is fixed automatically (the original fixes the first vertex of the
        /// first triangle). Components with no fixed vertex always get one added.
        /// </summary>
        public int[]? FixedVertices { get; set; }

        /// <summary>Relative residual tolerance of the linear solver. Default 1e-10.</summary>
        public double SolverTolerance { get; set; } = 1e-10;

        /// <summary>Maximum conjugate-gradient iterations per linear solve. Default 20000.</summary>
        public int SolverMaxIterations { get; set; } = 20000;
    }
}
