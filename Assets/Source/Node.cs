using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;

/// <summary>
/// Basic node class. Represents a 3D particle of the system.
/// </summary>
public class Node : MonoBehaviour
{
    #region InEditorVariables

    public bool Fixed;
    public float Mass;

    #endregion

    public int index;
    public Vector3 Pos;
    public Vector3 Pos0;
    public Vector3 Vel;

    private Simulable simul;

    private bool init = false;

    void Update()
    {
        if (init)
        {
            if (!Fixed)
            {
                updateUnityFromNode();
            }
            else
            {
                updateNodeFromUnity();
            }
        }
    }

    private void OnDrawGizmos()
    {
        if (this.Fixed)
        {
            Gizmos.color = Color.red;
            Gizmos.DrawSphere(Pos, this.simul.GetScale()*0.01f);
        }
        else
        {
            Gizmos.color = Color.green;
            Gizmos.DrawSphere(Pos, this.simul.GetScale()*0.01f);
        }
    }

    /// <summary>
    /// Sets unity transform position to simulation position.
    /// </summary>
    public void updateUnityFromNode()
    {
        transform.position = Pos;
    }

    /// <summary>
    /// Sets simulation position to unity transform position.
    /// </summary>
    public void updateNodeFromUnity()
    {
        Pos = transform.position;
    }

    /// <summary>
    /// Stablishes current position as the rest configuration.
    /// </summary>
    public void setCurrentAsPosRest()
    {
        Pos0 = Pos;
    }

    /// <summary>
    /// Initializes the simulation of the particle.
    /// </summary>
    /// <param name="idx">The global index from which it is assembled</param>
    /// <param name="sim">The simulable object that owns this particle</param>
    public void Initialize(int idx, Simulable sim)
    {
        index = idx;
        simul = sim;

        this.Mass = 0.0f;

        updateNodeFromUnity();
        setCurrentAsPosRest();

        this.Vel = new Vector3(0.0f, 0.0f, 0.0f);

        this.init = true;
    }

    /// <summary>
    /// Computes and returns the potential energy produced by this particle.
    /// Usually this is restricted to the work done by conservative external
    /// nodal forces (e.g., gravity).
    /// </summary>
    /// <returns>The nodal potential energy</returns>
    public double GetEnergy()
    {
        return -Vector3.Dot(Mass * simul.Gravity, Pos);
    }

    /// <summary>
    /// Computes and assembles the forces produced by this particle. Usually
    /// this is restricted to conservative external nodal forces (e.g., gravity),
    /// and damping forces (drag).
    /// </summary>
    /// <param name="force">The globlal force vector where to assemble local force</param>
    public void AddForce(VectorXD force)
    {
        Vector3 nodeForce = Vector3.zero;
        nodeForce += Mass * simul.Gravity; // Gravity force
        nodeForce += simul.Drag * Mass*Vel; // Drag force
        force.SetSubVector(index, 3, Utils.ToVectorXD(nodeForce) + force.SubVector(index, 3));
    }

    /// <summary>
    /// Computes and assembles the Jacobians of the forces produced by this particle.
    /// Usually this is restricted to conservative external nodal forces (e.g., gravity)
    /// and damping forces (drag).
    /// </summary>
    /// <param name="DfDx">The global Jacobian matrix where to assemble local DfDx</param>
    /// <param name="DfDv">The global Jacobian matrix where to assemble local DfDv</param>
    public void AddJacobian(MatrixXD DfDx, MatrixXD DfDv)
    {
        DfDv.SetSubMatrix(index, index, -simul.Drag * Mass * DenseMatrixXD.CreateIdentity(3) + DfDv.SubMatrix(index, 3, index, 3));
    }

    /// <summary>
    /// Computes and assembles the Jacobians of the forces produced by this particle.
    /// Usually this is restricted to conservative external nodal forces (e.g., gravity)
    /// and damping forces (drag).
    /// </summary>
    /// <param name="DfDx">The global Jacobian matrix where to assemble local DfDx</param>
    /// <param name="DfDv">The global Jacobian matrix where to assemble local DfDv</param>
    public void AddJacobian(List<Triplet> DfDx, List<Triplet> DfDv)
    {
        MatrixXD localDfDv = -simul.Drag * Mass * DenseMatrixXD.CreateIdentity(3);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                DfDv.Add(new Triplet(index + i, index + j, localDfDv[i,j]));
    }

    /// <summary>
    /// Assembles the position of the particle into the global position vector.
    /// </summary>
    /// <param name="pos">The global position vector where to assemble local position</param>
    public void GetPosition(VectorXD pos)
    {
        pos[index + 0] = Pos.x;
        pos[index + 1] = Pos.y;
        pos[index + 2] = Pos.z;
    }

    /// <summary>
    /// Assembles the velocity of the particle into the global velocity vector.
    /// </summary>
    /// <param name="pos">The global velocity vector where to assemble local velocity</param>
    public void GetVelocity(VectorXD vel)
    {
        vel[index + 0] = Vel.x;
        vel[index + 1] = Vel.y;
        vel[index + 2] = Vel.z;
    }

    /// <summary>
    /// Sets the position of the particle from the global position vector.
    /// </summary>
    /// <param name="pos">The global position vector from where to get the position of the particle</param>
    public void SetPosition(VectorXD pos)
    {
        Pos = new Vector3((float)pos[index + 0], (float)pos[index + 1], (float)pos[index + 2]);
    }

    /// <summary>
    /// Sets the velocity of the particle from the global velocity vector.
    /// </summary>
    /// <param name="pos">The global velocity vector from where to get the velocity of the particle</param>
    public void SetVelocity(VectorXD vel)
    {
        Vel = new Vector3((float)vel[index + 0], (float)vel[index + 1], (float)vel[index + 2]);
    }

    /// <summary>
    /// Assembles the mass of the particle into the global mass matrix.
    /// </summary>
    /// <param name="mass">The global mass matrix where to assemble the mass of the particle</param>
    public void GetMass(MatrixXD mass)
    {
        mass[index + 0, index + 0] = Mass;
        mass[index + 1, index + 1] = Mass;
        mass[index + 2, index + 2] = Mass;
    }

    /// <summary>
    /// Assembles the mass of the particle into the global mass matrix.
    /// </summary>
    /// <param name="mass">The global mass matrix where to assemble the mass of the particle</param>
    public void GetMass(List<Triplet> mass)
    {
        mass.Add(new Triplet(index + 0, index + 0, Mass));
        mass.Add(new Triplet(index + 1, index + 1, Mass));
        mass.Add(new Triplet(index + 2, index + 2, Mass));
    }

    /// <summary>
    /// Assembles the mass inverse of the particle into the global mass matrix.
    /// </summary>
    /// <param name="mass">The global mass inverse matrix where to assemble the mass of the particle</param>
    public void GetMassInverse(MatrixXD massInv)
    {
        massInv[index + 0, index + 0] = 1.0 / Mass;
        massInv[index + 1, index + 1] = 1.0 / Mass;
        massInv[index + 2, index + 2] = 1.0 / Mass;
    }

    /// <summary>
    /// Assembles the mass of the particle into the global mass matrix.
    /// </summary>
    /// <param name="mass">The global mass matrix where to assemble the mass of the particle</param>
    public void GetMassInverse(List<Triplet> massInv)
    {
        massInv.Add(new Triplet(index + 0, index + 0, 1.0/Mass));
        massInv.Add(new Triplet(index + 1, index + 1, 1.0/Mass));
        massInv.Add(new Triplet(index + 2, index + 2, 1.0/Mass));
    }

    /// <summary>
    /// If the particle is fixed, sets the positions of the specified
    /// global vector corresponding to the particle to true, in order
    /// to fix those degrees-of-freedom during simulation.
    /// </summary>
    public void SetFixedStencil(bool[] stencil)
    {
        if (Fixed)
        {
            stencil[index + 0] = true;
            stencil[index + 1] = true;
            stencil[index + 2] = true;
        }
    }

    /// <summary>
    /// If the particle is fixed, sets the positions of the specified
    /// global vector corresponding to the particle to 0, in order to
    /// fix those particles during simulation.
    /// </summary>
    /// <param name="v">The global vector to fix</param>
    public void FixVector(VectorXD v)
    {
        if (Fixed)
        {
            v[index + 0] = 0.0;
            v[index + 1] = 0.0;
            v[index + 2] = 0.0;
        }
    }

    /// <summary>
    /// If the particle is fixed, sets the rows/cols of the specified
    /// global matrix corresponding to the particle to 0, in order to
    /// fix those particles during simulation.
    /// </summary>
    /// <param name="M">The global matrix to fix</param>
    public void FixMatrix(MatrixXD M)
    {
        if (Fixed)
        {
            for (int i = 0; i < M.RowCount; i++)
            {
                M[index + 0, i] = 0.0;
                M[index + 1, i] = 0.0;
                M[index + 2, i] = 0.0;
                M[i, index + 0] = 0.0;
                M[i, index + 1] = 0.0;
                M[i, index + 2] = 0.0;
            }
            M[index + 0, index + 0] = 1.0;
            M[index + 1, index + 1] = 1.0;
            M[index + 2, index + 2] = 1.0;
        }
    }

}
