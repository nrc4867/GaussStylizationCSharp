// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

using GaussStylization.Internal;

namespace GaussStylizationCSharp.Tests;

public class MeshTests
{
    [Fact]
    public void Topology_ClosedMeshCounts()
    {
        var (v, f) = TestMeshes.IcoSphere(1); // 42 vertices, 80 faces, 120 edges
        var topo = new MeshTopology(v, f);
        Assert.Equal(42, topo.VertexCount);
        Assert.Equal(80, topo.FaceCount);
        Assert.Equal(120, topo.EdgeCount);
        for (int e = 0; e < topo.EdgeCount; e++)
        {
            Assert.NotEqual(-1, topo.EdgeFaces[2 * e]);
            Assert.NotEqual(-1, topo.EdgeFaces[2 * e + 1]);
        }
    }

    [Fact]
    public void Topology_GridHasBoundaryEdges()
    {
        var (v, f) = TestMeshes.Grid(4); // 3x3 quads
        var topo = new MeshTopology(v, f);
        int boundary = 0;
        for (int e = 0; e < topo.EdgeCount; e++) if (topo.EdgeFaces[2 * e + 1] == -1) boundary++;
        Assert.Equal(12, boundary);
    }

    [Fact]
    public void Topology_FaceEdgesAreOppositeTheVertex()
    {
        var (v, f) = TestMeshes.IcoSphere(0);
        var topo = new MeshTopology(v, f);
        for (int face = 0; face < topo.FaceCount; face++)
            for (int j = 0; j < 3; j++)
            {
                int e = topo.FaceEdges[3 * face + j];
                int a = topo.EdgeVertices[2 * e], b = topo.EdgeVertices[2 * e + 1];
                int opposite = f[3 * face + j];
                Assert.NotEqual(opposite, a);
                Assert.NotEqual(opposite, b);
            }
    }

    [Fact]
    public void Topology_NonManifoldEdgeThrows()
    {
        double[] v = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1 };
        int[] f = { 0, 1, 2, 0, 1, 3, 0, 1, 4 };
        Assert.Throws<ArgumentException>(() => new MeshTopology(v, f));
    }

    [Fact]
    public void CotangentLaplacian_AnnihilatesLinearFunctionsOnAFlatGrid()
    {
        int n = 6;
        var (v, f) = TestMeshes.Grid(n);
        var topo = new MeshTopology(v, f);
        var fn = new double[n * n];
        for (int i = 0; i < n * n; i++) fn[i] = 3 * v[3 * i] - 2 * v[3 * i + 1] + 1;
        var lap = new double[n * n];
        for (int e = 0; e < topo.EdgeCount; e++)
        {
            int a = topo.EdgeVertices[2 * e], b = topo.EdgeVertices[2 * e + 1];
            double w = topo.EdgeWeight[e];
            lap[a] += w * (fn[b] - fn[a]);
            lap[b] += w * (fn[a] - fn[b]);
        }
        for (int y = 1; y < n - 1; y++)
            for (int x = 1; x < n - 1; x++)
                Assert.Equal(0.0, lap[y * n + x], 10);
    }

    [Fact]
    public void CotangentWeights_AreOneHalfForRightIsoscelesGridEdges()
    {
        var (v, f) = TestMeshes.Grid(3);
        var topo = new MeshTopology(v, f);
        // The diagonal of each quad has weight 0; axis-aligned interior edges have weight 1.
        int zero = 0, one = 0;
        for (int e = 0; e < topo.EdgeCount; e++)
        {
            double w = topo.EdgeWeight[e];
            if (Math.Abs(w) < 1e-12) zero++;
            else if (Math.Abs(w - 1) < 1e-12) one++;
        }
        Assert.Equal(4, zero);
        Assert.True(one >= 4);
    }

    [Fact]
    public void FixedVertices_PinOneVertexPerComponentAndUnusedVertices()
    {
        // two disjoint triangles plus an unreferenced vertex
        double[] v = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 5, 0, 0, 6, 0, 0, 5, 1, 0, 9, 9, 9 };
        int[] f = { 0, 1, 2, 3, 4, 5 };
        var topo = new MeshTopology(v, f);
        bool[] pinned = topo.ComputeFixedVertices(null);
        Assert.True(pinned[0]);
        Assert.True(pinned[3]);
        Assert.True(pinned[6]);
        Assert.Equal(3, pinned.Count(p => p));
    }

    [Fact]
    public void LaplacianSolver_ReproducesHarmonicFunctionFromBoundaryValues()
    {
        int n = 7;
        var (v, f) = TestMeshes.Grid(n);
        var topo = new MeshTopology(v, f);
        var isFixed = new bool[n * n];
        var values = new double[n * n];
        for (int y = 0; y < n; y++)
            for (int x = 0; x < n; x++)
            {
                int i = y * n + x;
                values[i] = 2 * x - 3 * y + 5;
                isFixed[i] = x == 0 || y == 0 || x == n - 1 || y == n - 1;
            }
        var solver = new LaplacianSolver(topo, isFixed, 1e-12, 5000);
        var sol = new double[n * n]; // zero initial guess
        solver.Solve(new double[n * n], values, sol);
        Assert.True(TestMeshes.MaxAbsDiff(sol, values) < 1e-8);
    }
}
