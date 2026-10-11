// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

using System;

namespace GaussStylization.Internal
{
    internal readonly struct Vec3
    {
        public readonly double X, Y, Z;

        public Vec3(double x, double y, double z) { X = x; Y = y; Z = z; }

        public static Vec3 operator +(Vec3 a, Vec3 b) => new Vec3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        public static Vec3 operator -(Vec3 a, Vec3 b) => new Vec3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        public static Vec3 operator -(Vec3 a) => new Vec3(-a.X, -a.Y, -a.Z);
        public static Vec3 operator *(double s, Vec3 a) => new Vec3(s * a.X, s * a.Y, s * a.Z);
        public static Vec3 operator *(Vec3 a, double s) => new Vec3(s * a.X, s * a.Y, s * a.Z);

        public double Dot(Vec3 b) => X * b.X + Y * b.Y + Z * b.Z;
        public Vec3 Cross(Vec3 b) => new Vec3(Y * b.Z - Z * b.Y, Z * b.X - X * b.Z, X * b.Y - Y * b.X);
        public double Norm => Math.Sqrt(Dot(this));
        public double SquaredNorm => Dot(this);

        public Vec3 Normalized()
        {
            double n = Norm;
            return n > 0 ? this * (1.0 / n) : this;
        }
    }
}
