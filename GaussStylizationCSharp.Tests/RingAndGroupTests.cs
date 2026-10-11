// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

using GaussStylization;
using GaussStylization.Internal;

namespace GaussStylizationCSharp.Tests;

public class RingAndGroupTests
{
    private static readonly Vec3 Y = new Vec3(0, 1, 0);

    // ---- g function with rings -----------------------------------------------------------------

    [Fact]
    public void GFunction_RingOnly_IsOneOnTheRing()
    {
        var g = new GFunction(Array.Empty<Vec3>(), 4.0, new[] { 0.3 }, Y, 0.5);
        double s = Math.Sqrt(1 - 0.3 * 0.3);
        Assert.Equal(1.0, g.Value(new Vec3(s, 0.3, 0)), 9);
        Assert.Equal(1.0, g.Value(new Vec3(0, 0.3, -s)), 9);
        Assert.True(g.Value(new Vec3(1, 0, 0)) < 1.0);
    }

    [Fact]
    public void GFunction_RingValueMatchesTheFormula()
    {
        double sigma = 3.0, c = 0.5, d = 0.2;
        var n = new Vec3(0, 1, 0);
        var g = new GFunction(new[] { n }, sigma, new[] { d }, Y, c);
        double wn = 1.0 / Math.Exp(sigma);               // g(n)=1 for a single normal
        double wr = 1.0 / Math.Exp(sigma);               // ring weight: exp(sigma(1-0)) wr = 1
        var x = new Vec3(0.6, 0.64, 0.48).Normalized();
        double expected = c * wn * Math.Exp(sigma * n.Dot(x)) + wr * Math.Exp(sigma * (1 - Math.Pow(x.Y - d, 2)));
        Assert.Equal(expected, g.Value(x), 12);
    }

    [Fact]
    public void GFunction_WithRings_GradientAndHessianMatchFiniteDifferences()
    {
        var dirs = new[] { new Vec3(0, 1, 0), new Vec3(0, -1, 0), new Vec3(1, 0, 0) };
        var axis = new Vec3(1, 2, -1).Normalized();
        var g = new GFunction(dirs, 2.5, new[] { 0.0, 0.6 }, axis, 0.7);
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
            Assert.True(((gp - gm) * (1 / (2 * h)) - hess.Column(i)).Norm < 1e-5);
        }
    }

    [Fact]
    public void GFunction_DuplicateRingHeights_Throw()
    {
        Assert.Throws<ArgumentException>(() => new GFunction(Array.Empty<Vec3>(), 4.0, new[] { 0.2, 0.2 }, Y, 0.5));
    }

    // ---- rings through the public API ------------------------------------------------------------

    /// <summary>Mean of |n_a| (1 - |n_a|): zero when every normal is parallel or perpendicular to the axis.</summary>
    private static double Cylindricity(double[] v, int[] f, double[] axis)
    {
        double sum = 0;
        int nf = f.Length / 3;
        for (int i = 0; i < nf; i++)
        {
            Vec3 p0 = MeshTopology.Position(v, f[3 * i]);
            Vec3 n = (MeshTopology.Position(v, f[3 * i + 1]) - p0).Cross(MeshTopology.Position(v, f[3 * i + 2]) - p0).Normalized();
            double a = Math.Abs(n.Dot(new Vec3(axis[0], axis[1], axis[2])));
            sum += a * (1 - a);
        }
        return sum / nf;
    }

    [Fact]
    public void CylinderStyle_FlattensNormalsOntoCapsAndWalls()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        double before = Cylindricity(v, f, new double[] { 0, 1, 0 });
        var result = GaussStylizer.Stylize(v, f, new[] { GaussStyle.Cylinder() }, null, null, null, 40);
        double after = Cylindricity(result, f, new double[] { 0, 1, 0 });
        Assert.All(result, x => Assert.True(double.IsFinite(x)));
        Assert.True(after < before - 0.03, $"cylindricity {before:F4} -> {after:F4}");
    }

    [Fact]
    public void RingOnlyStyle_IsAccepted()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        var style = new GaussStyle { RingHeights = new[] { 0.0 } };
        var result = GaussStylizer.Stylize(v, f, new[] { style }, null, null, null, 5);
        Assert.All(result, x => Assert.True(double.IsFinite(x)));
        Assert.True(TestMeshes.MaxAbsDiff(v, result) > 1e-4);
    }

    [Fact]
    public void ConeStyle_Runs()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        var result = GaussStylizer.Stylize(v, f, new[] { GaussStyle.Cone() }, null, null, null, 20);
        Assert.All(result, x => Assert.True(double.IsFinite(x)));
        Assert.True(TestMeshes.MaxAbsDiff(v, result) > 1e-3);
    }

    [Fact]
    public void RingAxis_FollowsTheGridAxes()
    {
        // Rotating mesh + grid together must equal rotating the result: the ring axis is the grid's Y.
        var (v, f) = TestMeshes.IcoSphere(2);
        double a = 0.8, c = Math.Cos(a), s = Math.Sin(a);
        double[] r = { c, -s, 0, s, c, 0, 0, 0, 1 }; // rotation about Z
        double[] axes = { r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8] };
        // Cube directions keep the Newton systems full rank, so the comparison is not decided by roundoff ties.
        var styles = new[] { new GaussStyle(GaussStylizer.CubeDirections()) { RingHeights = new[] { 0.0, 0.6 } } };

        var viaGrid = GaussStylizer.Stylize(Apply(r, v), f, styles, null, null, axes, 8);
        var viaWorld = Apply(r, GaussStylizer.Stylize(v, f, styles, null, null, null, 8));
        Assert.True(TestMeshes.MaxAbsDiff(viaGrid, viaWorld) < 1e-5, $"equivariance diff {TestMeshes.MaxAbsDiff(viaGrid, viaWorld)}");

        // and a cylinder along a different axis really is a different result
        var other = GaussStylizer.Stylize(Apply(r, v), f, styles, null, null, null, 8);
        Assert.True(TestMeshes.MaxAbsDiff(viaGrid, other) > 1e-3, $"other diff {TestMeshes.MaxAbsDiff(viaGrid, other)}");
    }

    // ---- style groups ---------------------------------------------------------------------------

    [Fact]
    public void SingleGroup_MatchesTheSingleStyleOverload()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        var simple = GaussStylizer.Stylize(v, f, GaussStylizer.CubeDirections(), null, null, 6);
        var groups = GaussStylizer.Stylize(v, f, new[] { GaussStyle.Cube() }, new int[f.Length / 3], null, null, 6);
        Assert.Equal(simple, groups);
    }

    [Fact]
    public void FaceGroups_DriveEachRegionTowardItsOwnStyle()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        int nf = f.Length / 3;
        double[] alongX = { 1, 0, 0, -1, 0, 0 };
        double[] alongZ = { 0, 0, 1, 0, 0, -1 };
        var styles = new[] { new GaussStyle(alongX), new GaussStyle(alongZ) };

        double CentroidY(int face) => (v[3 * f[3 * face] + 1] + v[3 * f[3 * face + 1] + 1] + v[3 * f[3 * face + 2] + 1]) / 3;
        var groups = new int[nf];
        for (int i = 0; i < nf; i++) groups[i] = CentroidY(i) > 0 ? 0 : 1;

        var result = GaussStylizer.Stylize(v, f, styles, groups, null, null, 40);

        bool Upper(int face) => groups[face] == 0 && CentroidY(face) > 0.3;
        bool Lower(int face) => groups[face] == 1 && CentroidY(face) < -0.3;
        double Gain(double[] dirs, Func<int, bool> region) =>
            TestMeshes.Alignment(result, f, dirs, region) - TestMeshes.Alignment(v, f, dirs, region);

        // each region improves toward its own directions much more than toward the other group's
        Assert.True(Gain(alongX, Upper) > 0.03, $"upper->X {Gain(alongX, Upper):F4}");
        Assert.True(Gain(alongX, Upper) > Gain(alongZ, Upper) + 0.02, $"upper X {Gain(alongX, Upper):F4} vs Z {Gain(alongZ, Upper):F4}");
        Assert.True(Gain(alongZ, Lower) > 0.03, $"lower->Z {Gain(alongZ, Lower):F4}");
        Assert.True(Gain(alongZ, Lower) > Gain(alongX, Lower) + 0.02, $"lower Z {Gain(alongZ, Lower):F4} vs X {Gain(alongX, Lower):F4}");
    }

    [Fact]
    public void GroupWithZeroMu_DoesNotMove()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        var result = GaussStylizer.Stylize(v, f, new[] { new GaussStyle(GaussStylizer.CubeDirections()) { Mu = 0 } }, null, null, null, 5);
        Assert.True(TestMeshes.MaxAbsDiff(v, result) < 1e-7);
    }

    [Fact]
    public void GroupsCanMixDiscreteAndRingStyles()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        int nf = f.Length / 3;
        var groups = new int[nf];
        for (int i = 0; i < nf; i++) groups[i] = i % 2;
        var result = GaussStylizer.Stylize(v, f, new[] { GaussStyle.Cube(), GaussStyle.Cylinder() }, groups, null, null, 5,
            new GaussStylizationOptions { InnerIterations = 2 });
        Assert.All(result, x => Assert.True(double.IsFinite(x)));
    }

    // ---- validation -----------------------------------------------------------------------------

    public static IEnumerable<object[]> InvalidInputs()
    {
        var (v, f) = TestMeshes.IcoSphere(0);
        int nf = f.Length / 3;
        GaussStyle Ok() => GaussStyle.Cube();
        IReadOnlyList<GaussStyle> One(GaussStyle s) => new[] { s };

        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Array.Empty<GaussStyle>(), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, new GaussStyle[] { null! }, null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle()), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle(new double[] { 1, 0 })), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle { RingHeights = new[] { 1.5 } }), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle { RingHeights = new[] { 0.1, 0.1 } }), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle(GaussStylizer.CubeDirections()) { Sigma = 0 }), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle(GaussStylizer.CubeDirections()) { Mu = -1 }), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(new GaussStyle(GaussStylizer.CubeDirections()) { RingContribution = -1 }), null, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(Ok()), new int[nf - 1], null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(Ok()), Enumerable.Repeat(1, nf).ToArray(), null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, One(Ok()), Enumerable.Repeat(-1, nf).ToArray(), null, null, 1)) };
    }

    [Theory]
    [MemberData(nameof(InvalidInputs))]
    public void InvalidStyleInput_ThrowsArgumentException(Func<double[]> call) => Assert.ThrowsAny<ArgumentException>(call);

    [Fact]
    public void NullStyles_ThrowArgumentNullException()
    {
        var (v, f) = TestMeshes.IcoSphere(0);
        Assert.Throws<ArgumentNullException>(() => GaussStylizer.Stylize(v, f, (IReadOnlyList<GaussStyle>)null!, null, null, null, 1));
    }

    private static double[] Apply(double[] r, double[] v)
    {
        var o = new double[v.Length];
        for (int i = 0; i < v.Length / 3; i++)
            for (int d = 0; d < 3; d++)
                o[3 * i + d] = r[3 * d] * v[3 * i] + r[3 * d + 1] * v[3 * i + 1] + r[3 * d + 2] * v[3 * i + 2];
        return o;
    }
}
