using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;

/// <summary>
/// Basic interface for any simulation model.
/// </summary>
public interface ISimulable
{
    /// <summary>
    /// Initializes the simulable.
    /// </summary>
    void Initialize(int i, PhysicsManager m, List<Fixer> fixers, List<GameObject> obstacles);

    /// <summary>
    /// Returns the number of model DOF.
    /// </summary>
    int GetNumDoFs();

    /// <summary>
    /// Writes position values into the position vector.
    /// </summary>
    void GetPosition(VectorXD position);

    /// <summary>
    /// Sets position values from the position vector.
    /// </summary>
    void SetPosition(VectorXD position);

    /// <summary>
    /// Writes velocity values into the velocity vector.
    /// </summary>
    void GetVelocity(VectorXD velocity);

    /// <summary>
    /// Sets velocity values from the velocity vector.
    /// </summary>
    void SetVelocity(VectorXD velocity);

    /// <summary>
    /// Gets total potential energy.
    /// </summary>
    double GetEnergy();

    /// <summary>
    /// Writes force values into the force vector.
    /// </summary>
    void GetForce(VectorXD force);

    /// <summary>
    /// Writes force jacobian values into the matrix.
    /// </summary>
    void GetJacobian(MatrixXD DfDx, MatrixXD DfDv);

    /// <summary>
    /// Writes force jacobian values into the matrix.
    /// </summary>
    void GetJacobian(List<Triplet> DfDx, List<Triplet> DfDv);

    /// <summary>
    /// Writes mass values into the mass matrix.
    /// </summary>
    void GetMass(MatrixXD mass);

    /// <summary>
    /// Writes mass values into the sparse mass matrix.
    /// </summary>
    void GetMass(List<Triplet> mass);

    /// <summary>
    /// Write inverse of mass values into the inverse mass matrix.
    /// </summary>
    void GetMassInverse(MatrixXD massInv);

    /// <summary>
    /// Write inverse of mass values into the inverse sparse mass matrix.
    /// </summary>
    void GetMassInverse(List<Triplet> mass);

    /// <summary>
    /// Modifies a boolean vector indicating whether or not
    /// each of the degrees-of-freedom of the simulable is
    /// fixed.
    /// </summary>
    void GetFixedStencil(bool[] stencil);

    /// <summary>
    /// Fix the vector considering the fixed dofs.
    /// This method overwrites with 0's those fixed 
    /// degrees of freedom.
    /// </summary>
    void FixVector(VectorXD v);

    /// <summary>
    /// Fix the matrix considering the fixed dofs.
    /// This method overwrites with 0's the columns and rows of the linear system at
    /// fixed dofs, and puts 1's in the diagonal. This is valid as long as the right 
    /// hand side of the LS is also zero at fixed DOF. 
    /// </summary>
    void FixMatrix(MatrixXD M);

    /// <summary>
    /// Returns the approximate scale of the simulable (for visualization).
    /// </summary>
    float GetScale();

}
