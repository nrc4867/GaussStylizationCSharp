// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

using System;

namespace GaussStylization
{
    /// <summary>
    /// One style: the set of surface normals a group of faces is attracted to. Corresponds to one
    /// "style group" of the original application (its <c>style_N</c>, <c>style_R</c>, <c>mu</c>,
    /// <c>sigma</c> and <c>caxiscontrib</c>).
    /// </summary>
    /// <remarks>
    /// All directions are expressed in the grid frame (see the <c>gridAxes</c> argument of
    /// <see cref="GaussStylizer"/>): the local Y axis is the axis the rings are built around.
    /// </remarks>
    public sealed class GaussStyle
    {
        /// <summary>Creates a style with no preferred directions and no rings; set at least one of them.</summary>
        public GaussStyle() { }

        /// <summary>Creates a style from discrete preferred directions.</summary>
        public GaussStyle(double[] directions) { Directions = directions; }

        /// <summary>
        /// Discrete preferred normals, flat xyz triples (any number, may be empty when rings are given).
        /// Normalised internally; must be distinct.
        /// </summary>
        public double[] Directions { get; set; } = Array.Empty<double>();

        /// <summary>
        /// Optional rings (the "semi-discrete" style of the original, used for cylinders and cones).
        /// Each value is a height <c>d</c> in [-1, 1] along the grid's Y axis; the ring is the set of
        /// unit normals <c>n</c> with <c>n.y = d</c>. A single ring at 0 allows every normal perpendicular to the axis
        /// (together with directions +-Y this describes a cylinder). Heights must be distinct.
        /// </summary>
        public double[]? RingHeights { get; set; }

        /// <summary>Strength of the style for the faces of this group. Must be &gt;= 0. Default 1.</summary>
        public double Mu { get; set; } = 1.0;

        /// <summary>
        /// Sharpness of the g function; larger values attract faces more tightly to the preferred
        /// directions. Must be &gt; 0. Default 4.
        /// </summary>
        public double Sigma { get; set; } = 4.0;

        /// <summary>
        /// Weight of the discrete directions relative to the rings (the original's <c>caxiscontrib</c>).
        /// Only used when <see cref="RingHeights"/> is not empty. Must be &gt;= 0. Default 0.5.
        /// </summary>
        public double RingContribution { get; set; } = 0.5;

        /// <summary>The six cube normals (+-X, +-Y, +-Z).</summary>
        public static GaussStyle Cube() => new GaussStyle(GaussStylizer.CubeDirections());

        /// <summary>Cylinder along the grid's Y axis: normals +-Y plus a ring of all normals perpendicular to Y.</summary>
        public static GaussStyle Cylinder() => new GaussStyle(new double[] { 0, 1, 0, 0, -1, 0 })
        {
            RingHeights = new[] { 0.0 },
        };

        /// <summary>Cone along the grid's Y axis, as in the original application: normal -Y plus a ring at height 0.6.</summary>
        public static GaussStyle Cone() => new GaussStyle(new double[] { 0, -1, 0 })
        {
            RingHeights = new[] { 0.6 },
        };
    }
}
