// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

using GaussStylization.Internal;

namespace GaussStylizationCSharp.Tests;

public class InternalMathTests
{
    private static Mat3 RotationZ(double a) => new Mat3
    {
        M00 = Math.Cos(a), M01 = -Math.Sin(a), M10 = Math.Sin(a), M11 = Math.Cos(a), M22 = 1,
    };

    private static Mat3 RotationX(double a) => new Mat3
    {
        M00 = 1, M11 = Math.Cos(a), M12 = -Math.Sin(a), M21 = Math.Sin(a), M22 = Math.Cos(a),
    };

    private static void AssertProperRotation(Mat3 r)
    {
        Assert.Equal(1.0, r.Determinant, 9);
        Mat3 should = r * r.Transposed();
        Assert.Equal(1.0, should.M00, 9); Assert.Equal(1.0, should.M11, 9); Assert.Equal(1.0, should.M22, 9);
        Assert.Equal(0.0, should.M01, 9); Assert.Equal(0.0, should.M02, 9); Assert.Equal(0.0, should.M12, 9);
    }

    [Fact]
    public void FitRotation_RecoversRotationFromScaledTranspose()
    {
        Mat3 q = RotationZ(0.7) * RotationX(-1.1);
        Mat3 r = Mat3.FitRotation(2.5 * q.Transposed()); // tr(R S) is maximal at R = Q
        Assert.True(Max(r - q) < 1e-9);
    }

    [Fact]
    public void FitRotation_ReturnsProperRotationForReflections()
    {
        var s = new Mat3 { M00 = 1, M11 = 2, M22 = -3 };
        AssertProperRotation(Mat3.FitRotation(s));
    }

    [Fact]
    public void FitRotation_HandlesRankDeficientInput()
    {
        var a = new Vec3(1, 2, 3);
        var b = new Vec3(-2, 0.5, 1);
        Mat3 r = Mat3.FitRotation(Mat3.Outer(a, b));
        AssertProperRotation(r);
        // R aligns b with a, so tr(R S) = |a||b|
        Mat3 rs = r * Mat3.Outer(a, b);
        Assert.Equal(a.Norm * b.Norm, rs.M00 + rs.M11 + rs.M22, 8);
    }

    [Fact]
    public void FitRotation_ZeroMatrixGivesIdentity()
    {
        Assert.True(Max(Mat3.FitRotation(default) - Mat3.Identity) == 0);
    }

    [Fact]
    public void SolveSymmetric_SolvesIndefiniteSystem()
    {
        var a = new Mat3 { M00 = 2, M01 = 1, M02 = 0, M10 = 1, M11 = -3, M12 = 0.5, M20 = 0, M21 = 0.5, M22 = 1 };
        var x = new Vec3(0.3, -1.2, 2.0);
        Vec3 sol = Mat3.SolveSymmetric(a, a * x);
        Assert.True((sol - x).Norm < 1e-12);
    }

    [Fact]
    public void SolveSymmetric_SingularSystemStaysFinite()
    {
        Mat3 a = Mat3.Outer(new Vec3(1, 1, 0), new Vec3(1, 1, 0));
        Vec3 sol = Mat3.SolveSymmetric(a, new Vec3(1, 1, 0));
        Assert.False(double.IsNaN(sol.X + sol.Y + sol.Z));
        Assert.True((a * sol - new Vec3(1, 1, 0)).Norm < 1e-9); // consistent system is still solved
    }

    [Fact]
    public void GFunction_IsOneAtEveryPreferredNormal()
    {
        var dirs = new[] { new Vec3(0, 0, 1), new Vec3(0, 1, 0), new Vec3(1, 0, 0), new Vec3(0, 0, -1) };
        var g = new GFunction(dirs, 4.0);
        foreach (var d in dirs) Assert.Equal(1.0, g.Value(d), 9);
    }

    [Fact]
    public void GFunction_GradientAndHessianMatchFiniteDifferences()
    {
        var dirs = new[] { new Vec3(0, 0, 1), new Vec3(0.6, 0.8, 0), new Vec3(-1, 0, 0) };
        var g = new GFunction(dirs, 3.0);
        var x = new Vec3(0.2, -0.5, 0.7);
        g.GradientAndHessian(x, out Vec3 grad, out Mat3 hess);

        const double h = 1e-6;
        for (int i = 0; i < 3; i++)
        {
            Vec3 e = i == 0 ? new Vec3(h, 0, 0) : i == 1 ? new Vec3(0, h, 0) : new Vec3(0, 0, h);
            double fd = (g.Value(x + e) - g.Value(x - e)) / (2 * h);
            Assert.Equal(fd, i == 0 ? grad.X : i == 1 ? grad.Y : grad.Z, 6);

            g.GradientAndHessian(x + e, out Vec3 gp, out _);
            g.GradientAndHessian(x - e, out Vec3 gm, out _);
            Vec3 col = (gp - gm) * (1 / (2 * h));
            Vec3 hcol = hess.Column(i);
            Assert.True((col - hcol).Norm < 1e-5);
        }
    }

    [Fact]
    public void GFunction_RejectsDuplicateDirections()
    {
        Assert.Throws<ArgumentException>(() => new GFunction(new[] { new Vec3(0, 0, 1), new Vec3(0, 0, 1) }, 4.0));
    }

    private static double Max(Mat3 m) => new[]
    {
        m.M00, m.M01, m.M02, m.M10, m.M11, m.M12, m.M20, m.M21, m.M22,
    }.Max(Math.Abs);
}
