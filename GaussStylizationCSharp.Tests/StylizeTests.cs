// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

using GaussStylization;

namespace GaussStylizationCSharp.Tests;

public class StylizeTests
{
    private static readonly double[] Cube = GaussStylizer.CubeDirections();

    [Fact]
    public void ZeroIterations_ReturnsUnchangedCopy()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        var original = (double[])v.Clone();
        var result = GaussStylizer.Stylize(v, f, Cube, null, null, 0);
        Assert.NotSame(v, result);
        Assert.Equal(original, result);
    }

    [Fact]
    public void InputArraysAreNotModified()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        var v0 = (double[])v.Clone();
        var f0 = (int[])f.Clone();
        GaussStylizer.Stylize(v, f, Cube, null, null, 3);
        Assert.Equal(v0, v);
        Assert.Equal(f0, f);
    }

    [Fact]
    public void ZeroWeights_LeaveTheMeshUnchanged()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        var weights = new double[v.Length / 3];
        var result = GaussStylizer.Stylize(v, f, Cube, weights, null, 5);
        Assert.True(TestMeshes.MaxAbsDiff(v, result) < 1e-7);
    }

    [Fact]
    public void ZeroLambda_LeavesTheMeshUnchanged()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        var result = GaussStylizer.Stylize(v, f, Cube, null, null, 5, new GaussStylizationOptions { Lambda = 0 });
        Assert.True(TestMeshes.MaxAbsDiff(v, result) < 1e-7);
    }

    [Fact]
    public void Stylizing_MovesFaceNormalsTowardThePreferredDirections()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        double before = TestMeshes.Alignment(v, f, Cube);
        var result = GaussStylizer.Stylize(v, f, Cube, null, null, 40);
        double after = TestMeshes.Alignment(result, f, Cube);

        Assert.All(result, x => Assert.True(double.IsFinite(x)));
        Assert.True(after > before + 0.03, $"alignment {before:F4} -> {after:F4}");
    }

    [Fact]
    public void MoreIterations_AlignAtLeastAsWell()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        double few = TestMeshes.Alignment(GaussStylizer.Stylize(v, f, Cube, null, null, 5), f, Cube);
        double many = TestMeshes.Alignment(GaussStylizer.Stylize(v, f, Cube, null, null, 40), f, Cube);
        Assert.True(many >= few - 1e-6, $"{few:F4} vs {many:F4}");
    }

    [Fact]
    public void FirstVertexOfFirstTriangle_StaysPinnedByDefault()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        var result = GaussStylizer.Stylize(v, f, Cube, null, null, 10);
        int pinned = f[0];
        for (int d = 0; d < 3; d++) Assert.Equal(v[3 * pinned + d], result[3 * pinned + d], 12);
    }

    [Fact]
    public void ExplicitFixedVertices_StayInPlace()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        int[] pins = { 3, 17, 30 };
        var result = GaussStylizer.Stylize(v, f, Cube, null, null, 10, new GaussStylizationOptions { FixedVertices = pins });
        foreach (int p in pins)
            for (int d = 0; d < 3; d++) Assert.Equal(v[3 * p + d], result[3 * p + d], 12);
    }

    [Fact]
    public void GridAxes_RotateThePreferredDirections()
    {
        // Stylizing a rotated mesh with grid axes R must equal rotating the result of stylizing the
        // original mesh with the world axes (the problem is rotation-equivariant about the origin).
        var (v, f) = TestMeshes.IcoSphere(2);
        double a = 0.6, c = Math.Cos(a), s = Math.Sin(a);
        // rotation about the (1,1,0)/sqrt2-ish axis built from two plane rotations
        double[] rz = { c, -s, 0, s, c, 0, 0, 0, 1 };
        double[] rx = { 1, 0, 0, 0, c, -s, 0, s, c };
        double[] r = Mul(rz, rx);
        // gridAxes rows are the axes = columns of R
        double[] axes = { r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8] };

        var rotatedV = Apply(r, v);
        var viaGrid = GaussStylizer.Stylize(rotatedV, f, Cube, null, axes, 8);
        var viaWorld = Apply(r, GaussStylizer.Stylize(v, f, Cube, null, null, 8));
        Assert.True(TestMeshes.MaxAbsDiff(viaGrid, viaWorld) < 1e-5, $"{TestMeshes.MaxAbsDiff(viaGrid, viaWorld)}");
    }

    [Fact]
    public void GridAxes_ChangeTheResult()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        double c = Math.Cos(0.5), s = Math.Sin(0.5);
        double[] axes = { c, s, 0, -s, c, 0, 0, 0, 1 };
        var plain = GaussStylizer.Stylize(v, f, Cube, null, null, 5);
        var turned = GaussStylizer.Stylize(v, f, Cube, null, axes, 5);
        Assert.True(TestMeshes.MaxAbsDiff(plain, turned) > 1e-3);
    }

    [Fact]
    public void VertexWeights_ConfineTheEffectToWeightedRegions()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        int nv = v.Length / 3;
        var weights = new double[nv];
        for (int i = 0; i < nv; i++) weights[i] = v[3 * i + 1] > 0 ? 1.0 : 0.0; // stylize the upper half only
        var result = GaussStylizer.Stylize(v, f, Cube, weights, null, 30);

        // Compare how much each region's face normals gathered toward the cube directions.
        bool Upper(int face) => (v[3 * f[3 * face] + 1] + v[3 * f[3 * face + 1] + 1] + v[3 * f[3 * face + 2] + 1]) / 3 > 0.3;
        bool Lower(int face) => (v[3 * f[3 * face] + 1] + v[3 * f[3 * face + 1] + 1] + v[3 * f[3 * face + 2] + 1]) / 3 < -0.3;
        double upperGain = TestMeshes.Alignment(result, f, Cube, Upper) - TestMeshes.Alignment(v, f, Cube, Upper);
        double lowerGain = TestMeshes.Alignment(result, f, Cube, Lower) - TestMeshes.Alignment(v, f, Cube, Lower);
        Assert.True(upperGain > 0.03, $"weighted region should align: gain {upperGain:F4}");
        Assert.True(upperGain > lowerGain + 0.02, $"upper gain {upperGain:F4} vs lower gain {lowerGain:F4}");
    }

    [Fact]
    public void Result_IsDeterministicOnLargeParallelMeshes()
    {
        var (v, f) = TestMeshes.IcoSphere(3); // 1280 faces: exercises the parallel code paths
        var a = GaussStylizer.Stylize(v, f, Cube, null, null, 3);
        var b = GaussStylizer.Stylize(v, f, Cube, null, null, 3);
        Assert.Equal(a, b);
    }

    [Fact]
    public void SinglePrecisionOverload_MatchesDoublePrecision()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        var vf = v.Select(x => (float)x).ToArray();
        var df = Cube.Select(x => (float)x).ToArray();
        var single = GaussStylizer.Stylize(vf, f, df, null, null, 5);
        var dbl = GaussStylizer.Stylize(vf.Select(x => (double)x).ToArray(), f, Cube, null, null, 5);
        for (int i = 0; i < single.Length; i++) Assert.Equal(dbl[i], single[i], 5);
    }

    [Fact]
    public void CustomDirections_AreSupported()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        double[] up = { 0, 1, 0, 0, -1, 0 }; // a "plane" style
        double before = TestMeshes.Alignment(v, f, up);
        var result = GaussStylizer.Stylize(v, f, up, null, null, 30, new GaussStylizationOptions { InnerIterations = 3 });
        Assert.True(TestMeshes.Alignment(result, f, up) > before);
    }

    [Fact]
    public void OpenMesh_WithBoundary_Works()
    {
        var (v, f) = TestMeshes.Grid(8);
        var result = GaussStylizer.Stylize(v, f, Cube, null, null, 5);
        Assert.All(result, x => Assert.True(double.IsFinite(x)));
    }

    [Fact]
    public void EmptyMesh_ReturnsEmpty()
    {
        Assert.Empty(GaussStylizer.Stylize(Array.Empty<double>(), Array.Empty<int>(), Cube, null, null, 3));
    }

    public static IEnumerable<object[]> InvalidInputs()
    {
        var (v, f) = TestMeshes.IcoSphere(0);
        int nv = v.Length / 3;
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(new double[4], f, Cube, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, new[] { 0, 1 }, Cube, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, new[] { 0, 1, 99 }, Cube, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, new[] { 0, 1, 1 }, Cube, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Array.Empty<double>(), null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, new double[] { 0, 0, 0 }, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, new double[] { 0, 0, 1, 0, 0, 2 }, null, null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, new double[nv - 1], null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, Enumerable.Repeat(-1.0, nv).ToArray(), null, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, null, new double[8], 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, null, new double[] { 2, 0, 0, 0, 1, 0, 0, 0, 1 }, 1)) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, null, null, 1, new GaussStylizationOptions { Lambda = -1 })) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, null, null, 1, new GaussStylizationOptions { InnerIterations = 0 })) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(v, f, Cube, null, null, 1, new GaussStylizationOptions { FixedVertices = new[] { nv } })) };
        yield return new object[] { new Func<double[]>(() => GaussStylizer.Stylize(new[] { double.NaN, 0, 0, 1, 0, 0, 0, 1, 0 }, new[] { 0, 1, 2 }, Cube, null, null, 1)) };
    }

    [Theory]
    [MemberData(nameof(InvalidInputs))]
    public void InvalidInput_ThrowsArgumentException(Func<double[]> call) => Assert.ThrowsAny<ArgumentException>(call);

    [Fact]
    public void NullArrays_ThrowArgumentNullException()
    {
        var (v, f) = TestMeshes.IcoSphere(0);
        Assert.Throws<ArgumentNullException>(() => GaussStylizer.Stylize(null!, f, Cube, null, null, 1));
        Assert.Throws<ArgumentNullException>(() => GaussStylizer.Stylize(v, null!, Cube, null, null, 1));
        Assert.Throws<ArgumentNullException>(() => GaussStylizer.Stylize(v, f, null!, null, null, 1));
    }

    [Fact]
    public void NegativeIterations_Throw()
    {
        var (v, f) = TestMeshes.IcoSphere(0);
        Assert.Throws<ArgumentOutOfRangeException>(() => GaussStylizer.Stylize(v, f, Cube, null, null, -1));
    }

    private static double Sq(double x) => x * x;

    private static double[] Mul(double[] a, double[] b)
    {
        var r = new double[9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++) r[3 * i + j] += a[3 * i + k] * b[3 * k + j];
        return r;
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
