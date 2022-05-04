using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;

/// <summary>
/// Basic linear/co-rotational FEM system.
/// </summary>
public class FEMSystem : Simulable
{
    /// <summary>
    /// Default constructor. All zero. 
    /// </summary>
    public FEMSystem()
    {
        Manager = null;
    }

    #region EditorVariables

    /// <summary>
    /// Whether the collisions are computed at simulated nodes
    /// or at embedded nodes. CAUTION: Selecting embedded nodes
    /// requieres to have an embedded mesh.
    /// </summary>
    public enum CollisionType
    {
        SimulatedNode,
        EmbeddedNode
    }

    public bool Corotational = false;  // Use co-rotational method
    public bool Plastic = false; // Use plasticity
    public float YieldTol = 0.2f; // Plasticity yield tolerance parameter
    public float YieldMax = 0.4f; // Plasticity yield maximum parameter
    public float Creep = 0.01f; // Plasticity creep parameters
    public float Young = 1.0e6f; // Young modulus
    public float Poisson = 0.35f; // Poisson's ratio
    public float AlphaDamping = 0.0f; // Alpha parameter for the Rayleight damping
    public float BetaDamping = 0.0f; // Beta parameter for the Rayleight damping
    public CollisionType CollType = CollisionType.SimulatedNode; // Collision type
    public float CollStiffness = 1000000; // The collision forces stiffness

    [HideInInspector] public List<GameObject> CollObstacles; // Collection of obstacles (planes)
    [HideInInspector] public List<Node> Nodes; // List of nodes of the system (set before initialize)
    [HideInInspector] public List<Tetrahedron> Tetras; // List of tetrahedrons of the system (set before initialize)

    #endregion

    private float scale = 1; // Approximate scale of the object to be used for visualization

    /// <summary>
    /// Stores the iso-parametric coordinate of an embedded point
    /// together with the current deformed position. The deformed
    /// position gets updated whenever the deformation of the
    /// tetrahedras changes.
    /// </summary>
    public class EmbeddedPoint
    {
        public VectorXD isoParam;
        public VectorXD position;
    }

    private Mesh embededMesh = null;
    private List<EmbeddedPoint> embeddedPoints;


    #region ISimulable

    /// <summary>
    /// Initializes the simulation of the FEM system.
    /// </summary>
    /// <param name="idx">The global index from which it is assembled</param>
    /// <param name="m">The PhysicsManager object that is simulating it</param>
    override public void Initialize(int idx, PhysicsManager m, List<Fixer> fixers, List<GameObject> obstacles)
    {
        base.Initialize(idx, m, fixers, obstacles);

        // Initialize nodes and tetrahedra

        float minX = Mathf.Infinity, minY = Mathf.Infinity, minZ = Mathf.Infinity;
        float maxX = -Mathf.Infinity, maxY = -Mathf.Infinity, maxZ = -Mathf.Infinity;

        for (int i = 0; i < Nodes.Count; ++i)
        {
            Nodes[i].Initialize(idx + 3 * i, this);
            for (int j = 0; j < fixers.Count; ++j)
                if (fixers[j].IsInside(Nodes[i].Pos))
                    Nodes[i].Fixed = true;

            if (Nodes[i].Pos.x < minX) minX = Nodes[i].Pos.x;
            if (Nodes[i].Pos.y < minY) minY = Nodes[i].Pos.y;
            if (Nodes[i].Pos.z < minZ) minZ = Nodes[i].Pos.z;
            if (Nodes[i].Pos.x > maxX) maxX = Nodes[i].Pos.x;
            if (Nodes[i].Pos.y > maxY) maxY = Nodes[i].Pos.y;
            if (Nodes[i].Pos.z > maxZ) maxZ = Nodes[i].Pos.z;
        }

        // Compute the approximate scale of the object

        this.scale = Mathf.Max(maxX - minX, maxY - minY, maxZ - minZ);

        // Intiliaze the tetrahedras of the system

        for (int i = 0; i < Tetras.Count; ++i)
            Tetras[i].Initialize(i, this);

        // If possible, embed the mesh

        if (this.embededMesh != null)
        {
            this.embeddedPoints = new List<EmbeddedPoint>();

            if (!this.EmbedMesh(this.embededMesh, this.embeddedPoints))
            {
                this.embeddedPoints = null;
                this.embededMesh = null;
            }
        }

        this.CollObstacles = obstacles;
    }

    public override int GetNumDoFs()
    {
        return 3 * Nodes.Count;
    }

    public override void GetPosition(VectorXD position)
    {
        foreach (Node node in Nodes)
        {
            node.GetPosition(position);
        }
    }

    public override void SetPosition(VectorXD position)
    {
        foreach (Node node in Nodes)
        {
            node.SetPosition(position);
        }

        if (Corotational)
        {
            foreach (Tetrahedron tetra in Tetras)
            {
                tetra.UpdateRotation();
            }
        }

        if (Plastic)
        {
            foreach (Tetrahedron tetra in Tetras)
            {
                tetra.UpdatePlasticity();
            }
        }

        if (this.embededMesh != null)
        {
            foreach (Tetrahedron tetra in Tetras)
            {
                tetra.UpdateEmbedded();
            }

            Vector3[] newVertices = new Vector3[this.embededMesh.vertexCount];
            for (int i = 0; i < this.embededMesh.vertexCount; ++i)
                newVertices[i] = Utils.ToVector3(this.embeddedPoints[i].position);

            this.embededMesh.vertices = newVertices;
            this.embededMesh.RecalculateBounds();
            this.embededMesh.RecalculateNormals();
            this.embededMesh.RecalculateTangents();
        }
    }

    public override void GetVelocity(VectorXD velocity)
    {
        foreach (Node node in Nodes)
        {
            node.GetVelocity(velocity);
        }
    }

    public override void SetVelocity(VectorXD velocity)
    {
        foreach (Node node in Nodes)
        {
            node.SetVelocity(velocity);
        }
    }

    public override double GetEnergy()
    {
        double energy = 0.0f;

        foreach (Node node in Nodes)
        {
            energy += node.GetEnergy();
        }

        foreach (Tetrahedron tetra in Tetras)
        {
            energy += tetra.GetEnergy();
        }

        return energy;
    }

    public override void GetForce(VectorXD force)
    {
        foreach (Node node in Nodes)
        {
            node.AddForce(force);
        }

        foreach (Tetrahedron tetra in Tetras)
        {
            tetra.GetForce(force);
        }

        this.AddCollisionForce(force);
    }

    public override void GetJacobian(MatrixXD dFdx, MatrixXD dFdv)
    {
        foreach (Node node in Nodes)
        {
            node.AddJacobian(dFdx, dFdv);
        }
        foreach (Tetrahedron tetra in Tetras)
        {
            tetra.GetJacobian(dFdx, dFdv);
        }

        this.AddCollisionJacobian(dFdx, dFdv);
    }

    public override void GetJacobian(List<Triplet> dFdx, List<Triplet> dFdv)
    {
        foreach (Node node in Nodes)
        {
            node.AddJacobian(dFdx, dFdv);
        }
        foreach (Tetrahedron tetra in Tetras)
        {
            tetra.GetJacobian(dFdx, dFdv);
        }

        this.AddCollisionJacobian(dFdx, dFdv);
    }

    public override void GetMass(MatrixXD mass)
    {
        foreach (Node node in Nodes)
        {
            node.GetMass(mass);
        }
    }

    public override void GetMass(List<Triplet> mass)
    {
        foreach (Node node in Nodes)
        {
            node.GetMass(mass);
        }
    }

    public override void GetMassInverse(MatrixXD massInv)
    {
        foreach (Node node in Nodes)
        {
            node.GetMassInverse(massInv);
        }
    }

    public override void GetMassInverse(List<Triplet> massInv)
    {
        foreach (Node node in Nodes)
        {
            node.GetMassInverse(massInv);
        }
    }

    public override void GetFixedStencil(bool[] stencil)
    {
        foreach (Node node in Nodes)
        {
            node.SetFixedStencil(stencil);
        }
    }

    public override void FixVector(VectorXD v)
    {
        foreach (Node node in Nodes)
        {
            node.FixVector(v);
        }
    }

    public override void FixMatrix(MatrixXD M)
    {
        foreach (Node node in Nodes)
        {
            node.FixMatrix(M);
        }
    }

    public override float GetScale() { return scale; }

    #endregion

    #region Collisions

    void AddCollisionForce(VectorXD force)
    {
        if (this.CollObstacles.Count != 0)
        {
            if (this.CollType == CollisionType.SimulatedNode)
                this.AddCollisionForceSimulated(force);
            else this.AddCollisionForceEmbedded(force);
        }
    }

    void AddCollisionJacobian(MatrixXD dfdx, MatrixXD dfdv)
    {
        if (this.CollObstacles.Count != 0)
        {
            if (this.CollType == CollisionType.SimulatedNode)
                this.AddCollisionJacobianSimulated(dfdx, dfdv);
            else this.AddCollisionJacobianEmbedded(dfdx, dfdv);
        }
    }

    void AddCollisionJacobian(List<Triplet> dfdx, List<Triplet> dfdv)
    {
        if (this.CollObstacles.Count != 0)
        {
            if (this.CollType == CollisionType.SimulatedNode)
                this.AddCollisionJacobianSimulated(dfdx, dfdv);
            else this.AddCollisionJacobianEmbedded(dfdx, dfdv);
        }
    }

    #endregion

    #region SimulatedNodeCollisions

    /// <summary>
    /// For each of the potential obstacles, assemble the corresponding
    /// collision forces due to the collision of the simulated nodes.
    /// </summary>
    void AddCollisionForceSimulated(VectorXD force)
    {
        for (int o = 0; o < this.CollObstacles.Count; ++o)
        {
            Vector3 up = new Vector3(0.0f, 1.0f, 0.0f);
            Vector3 p = this.CollObstacles[o].transform.position;
            Vector3 n = this.CollObstacles[o].transform.rotation * up;
            for (int i = 0; i < Nodes.Count; ++i)
            {
                float projDist = Vector3.Dot((Nodes[i].Pos - p), n);
                float k = this.CollStiffness;
                if (projDist < 0.01)
                {
                    force.SetSubVector(Nodes[i].index, 3, force.SubVector(Nodes[i].index, 3) + Utils.ToVectorXD(-k * projDist * n));
                }
            }
        }
    }

    /// <summary>
    /// For each of the potential obstacles, assemble the corresponding
    /// collision Jacobians due to the collision of the simulated nodes.
    /// </summary>
    void AddCollisionJacobianSimulated(MatrixXD dfdx, MatrixXD dfdv)
    {
        for (int o = 0; o < this.CollObstacles.Count; ++o)
        {
            Vector3 up = new Vector3(0.0f, 1.0f, 0.0f);
            Vector3 p = this.CollObstacles[o].transform.position;
            Vector3 n = this.CollObstacles[o].transform.rotation * up;
            for (int i = 0; i < Nodes.Count; ++i)
            {
                float projDist = Vector3.Dot((Nodes[i].Pos - p), n);
                float k = this.CollStiffness;
                if (projDist < 0.01)
                {
                    MatrixXD nnT = Utils.ToVectorXD(n).OuterProduct(Utils.ToVectorXD(n));
                    dfdx.SetSubMatrix(Nodes[i].index, Nodes[i].index, dfdx.SubMatrix(Nodes[i].index, 3, Nodes[i].index, 3) + -k * nnT);
                }
            }
        }
    }

    /// <summary>
    /// For each of the potential obstacles, assemble the corresponding
    /// collision Jacobians due to the collision of the simulated nodes
    /// (sparse version).
    /// </summary>
    void AddCollisionJacobianSimulated(List<Triplet> dfdx, List<Triplet> dfdv)
    {
        for (int o = 0; o < this.CollObstacles.Count; ++o)
        {
            Vector3 up = new Vector3(0.0f, 1.0f, 0.0f);
            Vector3 p = this.CollObstacles[o].transform.position;
            Vector3 n = this.CollObstacles[o].transform.rotation * up;
            for (int i = 0; i < Nodes.Count; ++i)
            {
                float projDist = Vector3.Dot((Nodes[i].Pos - p), n);
                float k = this.CollStiffness;
                if (projDist < 0.01)
                {
                    MatrixXD localDfDx = -k * Utils.ToVectorXD(n).OuterProduct(Utils.ToVectorXD(n));
                    for (int ii = 0; ii < 3; ++ii)
                        for (int jj = 0; jj < 3; ++jj)
                            dfdx.Add(new Triplet(Nodes[i].index + ii, Nodes[i].index + jj, localDfDx[ii, jj]));
                }
            }
        }
    }

    #endregion

    #region Embedding

    /// <summary>
    /// Sets a visualization mesh embedded in the simulation.
    /// This is not actually embedded in the tetrahedra until
    /// the Initialize method is called.
    /// </summary>
    public void SetEmbeddedMesh(Mesh mesh)
    {
        this.embededMesh = mesh;
    }

    /// <summary>
    /// Embeds the the specified mesh in the FEM system to be used for visualization 
    /// and contact purposes. It returns true if all the nodes of the mesh have been 
    /// succesfully embedded in a tetrahedron. For each vertex of the mesh, this method
    /// checks if the vertex is within a tetrahedron in the rest configuration and stores
    /// the corresponding embedded point in the specified list.
    /// </summary>
    /// <param name="mesh">The mesh to embed</param>
    /// <param name="embeddedPoints">The list of resulting embedded points</param>
    /// <returns>Returns true if all the mesh vertices are within the volume</returns>
    bool EmbedMesh(Mesh mesh, List<EmbeddedPoint> embeddedPoints)
    {
        for (int i = 0; i < mesh.vertexCount; ++i)
        {
            bool found = false;

            // COMPLETAR Parte 2
            EmbeddedPoint point=new EmbeddedPoint();
            point.position =Utils.ToVectorXD( mesh.vertices[i]+transform.position);
            for(int j = 0; j <Tetras.Count; j++)
            {
                if (Tetras[j].TryAddEmbedded(point))
                {
                    embeddedPoints.Add(point);
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                Debug.Log("WARNING: point " + i + " inside the volume!");
            }
        }

        return true;
    }

    /// <summary>
    /// For each of the potential obstacles, assemble the corresponding
    /// collision forces due to the collision of the embedded mesh vertices.
    /// </summary>
    void AddCollisionForceEmbedded(VectorXD force)
    {
        for (int o = 0; o < this.CollObstacles.Count; ++o)
        {
            // COMPLETAR Parte 3
            if (CollObstacles[o].CompareTag("Sphere"))
            {
                VectorXD p = Utils.ToVectorXD(CollObstacles[o].transform.position);
                double r = CollObstacles[o].GetComponent<SphereCollider>().radius;
                for (int i = 0; i < Tetras.Count; i++)
                {
                    Tetras[i].AddCollisionForceEmbeddedSphere(force, CollStiffness, p, r);
                }
            }
            else
            {
                VectorXD p = Utils.ToVectorXD(CollObstacles[o].transform.position);
                VectorXD n = Utils.ToVectorXD(CollObstacles[o].transform.rotation * new Vector3(0, 1, 0));
                for (int i = 0; i < Tetras.Count; i++)
                {
                    Tetras[i].AddCollisionForceEmbedded(force, CollStiffness, p, n);
                }
            }
            
        }
    }

    /// <summary>
    /// For each of the potential obstacles, assemble the corresponding
    /// collision Jacobians due to the collision of the embedded mesh vertices.
    /// </summary>
    void AddCollisionJacobianEmbedded(MatrixXD dfdx, MatrixXD dfdv)
    {
        for (int o = 0; o < this.CollObstacles.Count; ++o)
        {
            // COMPLETAR Parte 3
            if (CollObstacles[o].CompareTag("Sphere"))
            {
                VectorXD p = Utils.ToVectorXD(CollObstacles[o].transform.position);
                double r = CollObstacles[o].GetComponent<SphereCollider>().radius;
                for (int i = 0; i < Tetras.Count; i++)
                {
                    Tetras[i].AddCollisionJacobianEmbeddedSphere(dfdx, dfdv, CollStiffness, p, r);
                }
            }
            else
            {
                VectorXD p = Utils.ToVectorXD(CollObstacles[o].transform.position);
                VectorXD n = Utils.ToVectorXD(CollObstacles[o].transform.rotation * new Vector3(0, 1, 0));
                for (int i = 0; i < Tetras.Count; i++)
                {
                    Tetras[i].AddCollisionJacobianEmbedded(dfdx, dfdv, CollStiffness, p, n);
                }
            }
            
        }
    }

    /// <summary>
    /// For each of the potential obstacles, assemble the corresponding
    /// collision Jacobians due to the collision of the embedded mesh vertices.
    /// (sparse version).
    /// </summary>
    void AddCollisionJacobianEmbedded(List<Triplet> dfdx, List<Triplet> dfdv)
    {
        for (int o = 0; o < this.CollObstacles.Count; ++o)
        {
            // OPCIONAL Parte 3
            if (CollObstacles[o].CompareTag("Sphere"))
            {
                VectorXD p = Utils.ToVectorXD(CollObstacles[o].transform.position);
                float r = CollObstacles[o].GetComponent<SphereCollider>().radius;
                for (int i = 0; i < Tetras.Count; i++)
                {
                    Tetras[i].AddCollisionJacobianEmbeddedSphere(dfdx, dfdv, CollStiffness, p, r);
                }
            }
            else
            {
                VectorXD p = Utils.ToVectorXD(CollObstacles[o].transform.position);
                VectorXD n = Utils.ToVectorXD(CollObstacles[o].transform.rotation * new Vector3(0, 1, 0));
                for (int i = 0; i < Tetras.Count; i++)
                {
                    Tetras[i].AddCollisionJacobianEmbedded(dfdx, dfdv, CollStiffness, p, n);
                }
            }
        }
    }

    #endregion

}