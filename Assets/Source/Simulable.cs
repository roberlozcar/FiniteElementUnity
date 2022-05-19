using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;

/// <summary>
/// 
/// </summary>
public abstract class Simulable : MonoBehaviour, ISimulable
{
    /// <summary>
    /// Default constructor. All zero. 
    /// </summary>
    public Simulable()
    {
        Manager = null;
    }

    #region EditorVariables

    public PhysicsManager Manager = null;
    public bool TestForce = false;
    public bool TestDfDx = false;
    public float Density = 1.0f;
    public float Drag = 0.0f;
    public Vector3 Gravity = new Vector3(0.0f, -9.8f, 0.0f);

    #endregion

    #region OtherVariables
    protected int index;
    #endregion

    #region MonoBehaviour

    public void Start()
    {

    }

    public void Update()
    {

    }

    #endregion

    #region ISimulable

    public virtual void Initialize(int idx, PhysicsManager m, List<Fixer> fixers, List<GameObject> obstacles)
    {
        Manager = m;
        index = idx;
    }

    public abstract int GetNumDoFs();
    public abstract void GetPosition(VectorXD position);
    public abstract void SetPosition(VectorXD position);
    public abstract void GetVelocity(VectorXD velocity);
    public abstract void SetVelocity(VectorXD velocity);
    public abstract double GetEnergy();
    public abstract void GetForce(VectorXD force);
    public abstract void GetJacobian(MatrixXD DfDx, MatrixXD DfDv);
    public abstract void GetJacobian(List<Triplet> DfDx, List<Triplet> DfDv);
    public abstract void GetMass(MatrixXD mass);
    public abstract void GetMass(List<Triplet> mass);
    public abstract void GetMassInverse(MatrixXD massInv);
    public abstract void GetMassInverse(List<Triplet> massInv);
    public abstract void GetFixedStencil(bool[] stencil);
    public abstract void FixVector(VectorXD v);
    public abstract void FixMatrix(MatrixXD M);
    public abstract float GetScale();

    #endregion
}
