// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.
// Edge topology and cotangent weights follow libigl's edge_topology / cotmatrix_entries.

using System;

namespace GaussStylization.Internal
{
    /// <summary>
    /// Edge/face connectivity plus cotangent weights of a triangle mesh. Corresponds to the
    /// EVList / FEList / EFList / L data of gaussStylization's <c>gauss_style_data</c>.
    /// </summary>
    internal sealed class MeshTopology
    {
        public readonly int VertexCount;
        public readonly int FaceCount;
        public readonly int EdgeCount;

        /// <summary>3 vertex indices per face.</summary>
        public readonly int[] Faces;
        /// <summary>2 vertex indices per edge (smaller index first).</summary>
        public readonly int[] EdgeVertices;
        /// <summary>2 faces per edge; -1 for a missing (boundary) face.</summary>
        public readonly int[] EdgeFaces;
        /// <summary>3 edge indices per face; entry j is the edge opposite the face's vertex j.</summary>
        public readonly int[] FaceEdges;
        /// <summary>Per face, per edge (opposite vertex j): half the cotangent of the angle at vertex j.</summary>
        public readonly double[] FaceCot;
        /// <summary>Per edge cotangent weight w_ij (the off-diagonal entry of the cotangent Laplacian).</summary>
        public readonly double[] EdgeWeight;
        /// <summary>CSR vertex -> incident faces.</summary>
        public readonly int[] VertexFaceStart;
        public readonly int[] VertexFaces;

        public MeshTopology(double[] vertices, int[] faces)
        {
            VertexCount = vertices.Length / 3;
            FaceCount = faces.Length / 3;
            Faces = faces;

            // ---- edges -------------------------------------------------------------------
            int nh = 3 * FaceCount;
            var keys = new long[nh];
            var order = new int[nh];
            for (int f = 0; f < FaceCount; f++)
            {
                for (int j = 0; j < 3; j++)
                {
                    int a = faces[3 * f + (j + 1) % 3], b = faces[3 * f + (j + 2) % 3];
                    int lo = Math.Min(a, b), hi = Math.Max(a, b);
                    keys[3 * f + j] = (long)lo * VertexCount + hi;
                    order[3 * f + j] = 3 * f + j;
                }
            }
            Array.Sort(keys, order);

            var edgeVerts = new int[nh * 2];
            var edgeFaces = new int[nh * 2];
            FaceEdges = new int[nh];
            int ne = 0;
            for (int i = 0; i < nh;)
            {
                int run = 1;
                while (i + run < nh && keys[i + run] == keys[i]) run++;
                if (run > 2)
                    throw new ArgumentException(
                        "The mesh is not edge-manifold: an edge is shared by more than two triangles.");
                int h0 = order[i];
                int f0 = h0 / 3, j0 = h0 % 3;
                edgeVerts[2 * ne] = (int)(keys[i] / VertexCount);
                edgeVerts[2 * ne + 1] = (int)(keys[i] % VertexCount);
                edgeFaces[2 * ne] = f0;
                FaceEdges[h0] = ne;
                if (run == 2)
                {
                    int h1 = order[i + 1];
                    edgeFaces[2 * ne + 1] = h1 / 3;
                    FaceEdges[h1] = ne;
                }
                else
                {
                    edgeFaces[2 * ne + 1] = -1;
                }
                ne++;
                i += run;
            }
            EdgeCount = ne;
            EdgeVertices = new int[2 * ne];
            EdgeFaces = new int[2 * ne];
            Array.Copy(edgeVerts, EdgeVertices, 2 * ne);
            Array.Copy(edgeFaces, EdgeFaces, 2 * ne);

            // ---- cotangent weights -------------------------------------------------------
            FaceCot = new double[nh];
            EdgeWeight = new double[ne];
            for (int f = 0; f < FaceCount; f++)
            {
                Vec3 p0 = Position(vertices, faces[3 * f]);
                Vec3 p1 = Position(vertices, faces[3 * f + 1]);
                Vec3 p2 = Position(vertices, faces[3 * f + 2]);
                double l0 = (p2 - p1).SquaredNorm; // edge opposite vertex 0
                double l1 = (p0 - p2).SquaredNorm; // opposite vertex 1
                double l2 = (p1 - p0).SquaredNorm; // opposite vertex 2
                double dblA = (p1 - p0).Cross(p2 - p0).Norm;
                if (!(dblA > 0)) continue; // degenerate triangle: contributes no weight
                FaceCot[3 * f + 0] = (l1 + l2 - l0) / dblA / 4.0;
                FaceCot[3 * f + 1] = (l2 + l0 - l1) / dblA / 4.0;
                FaceCot[3 * f + 2] = (l0 + l1 - l2) / dblA / 4.0;
                for (int j = 0; j < 3; j++) EdgeWeight[FaceEdges[3 * f + j]] += FaceCot[3 * f + j];
            }

            // ---- vertex -> faces ---------------------------------------------------------
            VertexFaceStart = new int[VertexCount + 1];
            for (int i = 0; i < nh; i++) VertexFaceStart[faces[i] + 1]++;
            for (int v = 0; v < VertexCount; v++) VertexFaceStart[v + 1] += VertexFaceStart[v];
            VertexFaces = new int[nh];
            var cursor = new int[VertexCount];
            for (int f = 0; f < FaceCount; f++)
                for (int j = 0; j < 3; j++)
                {
                    int v = faces[3 * f + j];
                    VertexFaces[VertexFaceStart[v] + cursor[v]++] = f;
                }
        }

        public static Vec3 Position(double[] v, int i) => new Vec3(v[3 * i], v[3 * i + 1], v[3 * i + 2]);

        /// <summary>Unit face normals (zero vector for degenerate faces), like igl::per_face_normals.</summary>
        public Vec3[] FaceNormals(double[] positions)
        {
            var n = new Vec3[FaceCount];
            for (int f = 0; f < FaceCount; f++)
            {
                Vec3 p0 = Position(positions, Faces[3 * f]);
                Vec3 c = (Position(positions, Faces[3 * f + 1]) - p0).Cross(Position(positions, Faces[3 * f + 2]) - p0);
                double r = c.Norm;
                n[f] = r == 0 ? new Vec3(0, 0, 0) : c * (1.0 / r);
            }
            return n;
        }

        /// <summary>
        /// Vertices that are held at their input position. The Laplacian system is singular without
        /// at least one fixed vertex per connected component, so one is added for every component that
        /// has none (the original code fixes only <c>F(0,0)</c>). Vertices used by no face are fixed too.
        /// </summary>
        public bool[] ComputeFixedVertices(int[]? userFixed)
        {
            var fixedV = new bool[VertexCount];
            if (userFixed != null) foreach (int v in userFixed) fixedV[v] = true;

            var used = new bool[VertexCount];
            var parent = new int[VertexCount];
            for (int i = 0; i < VertexCount; i++) parent[i] = i;
            int Find(int x)
            {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            }
            for (int f = 0; f < FaceCount; f++)
            {
                int a = Faces[3 * f], b = Faces[3 * f + 1], c = Faces[3 * f + 2];
                used[a] = used[b] = used[c] = true;
                parent[Find(b)] = Find(a);
                parent[Find(c)] = Find(a);
            }
            for (int v = 0; v < VertexCount; v++) if (!used[v]) fixedV[v] = true;

            var componentPinned = new bool[VertexCount];
            for (int v = 0; v < VertexCount; v++) if (used[v] && fixedV[v]) componentPinned[Find(v)] = true;
            for (int f = 0; f < FaceCount; f++)
            {
                int v = Faces[3 * f];
                int root = Find(v);
                if (!componentPinned[root]) { componentPinned[root] = true; fixedV[v] = true; }
            }
            return fixedV;
        }
    }
}
