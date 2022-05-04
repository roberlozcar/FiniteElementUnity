using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using MathNet.Numerics.LinearAlgebra.Factorization;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;

/// <summary>
/// Basic Tetrahedron class. Represents a single element of a linear
/// FEM discretization. Computes and stores the rest-state-depending
/// quantities and computes the corresponding strain energy, forces
/// and Jacobians necessary for simulation. It also allows to use
/// the well-known co-rotational model.
/// </summary>
public class Tetrahedron : MonoBehaviour
{
    #region InEditorVariables

    public List<Node> Nodes;

    #endregion

    private int TetIndex; // Index of this tetrahedron (not used)

    private FEMSystem Simul; // Simulable using this tetrahedron

    private MatrixXD Rest2Iso; // Basis-change matrix from rest to iso-parametric
    private MatrixXD Iso2Rest; // Basis-change matrix from iso-parametric to rest
    private MatrixXD FRestToMatrix; // 3x4 version of DxDX
    private MatrixXD FRestToVector; // 9x12 version of DxDX 
    private double Volume0; // Volume at the rest configuration

    private MatrixXD FPlastic; // 3x3 Plastic part of DxDX

    private MatrixXD RotationOne; // 3x3 current rotation
    private MatrixXD RotationAll; // 12x12 current rotations

    private double firstlame;
    private double secondlame;

    /// <summary>
    /// List of the points embedded in this tetrahedron
    /// </summary>
    private List<FEMSystem.EmbeddedPoint> embeddedPoints = new List<FEMSystem.EmbeddedPoint>();

    #region Initialization

    /// <summary>
    /// Initializes the simulation of the tetrahedron.
    /// </summary>
    /// <param name="idx">The global index from which it is assembled</param>
    /// <param name="sim">The simulable object that owns this particle</param>
    public void Initialize(int idx, FEMSystem sim)
    {
        TetIndex = idx;
        Simul = sim;

        this.SetCurrentAsRest();

        this.embeddedPoints.Clear();
    }

    /// <summary>
    /// Sets the current state as the rest state. Then, calls the
    /// necessary methods to compute and store the quantities that
    /// only depend on the rest configuration.
    /// </summary>
    public void SetCurrentAsRest()
    {
        // Set rest positions

        for (int i = 0; i < 4; ++i)
            Nodes[i].setCurrentAsPosRest();

        // Recompute transformations to iso-parametric space

        this.Iso2Rest = this.ComputeIsoparametric2Rest();
        this.Rest2Iso = this.ComputeRest2Isoparametric();

        // Precompute the part of the deformation gradient
        // that only depends on the undeformed configuration.

        this.PrecomputeFRest();

        // Set volume at rest

        this.Volume0 = this.ComputeVolume();

        // Compute and distribute mass

        float mass = (float)this.Volume0 * Simul.Density;
        for (int i = 0; i < 4; ++i)
            Nodes[i].Mass += mass / 4.0f;

        // Compute initial rotation

        this.RotationOne = DenseMatrixXD.CreateIdentity(3);
        this.RotationAll = DenseMatrixXD.CreateIdentity(12);

        // Initialize plastic deformation

        this.FPlastic = DenseMatrixXD.Create(3, 3, 0);

        firstlame=ComputeLameParameterFirst();
        secondlame=ComputeLameParameterSecond();
    }

    #endregion

    #region BasicMagnitudes

    /// <summary>
    /// Compute the volume of the tetrahedron in the current state.
    /// </summary>
    /// <returns>The volume of the tetrahedron</returns>
    public double ComputeVolume()
    {
        Vector3 e0 = this.Nodes[1].Pos - this.Nodes[0].Pos;
        Vector3 e1 = this.Nodes[2].Pos - this.Nodes[0].Pos;
        Vector3 e2 = this.Nodes[3].Pos - this.Nodes[0].Pos;
        MatrixXD mV = new DenseMatrixXD(3);
        mV.SetColumn(0, Utils.ToVectorXD(e0));
        mV.SetColumn(1, Utils.ToVectorXD(e1));
        mV.SetColumn(2, Utils.ToVectorXD(e2));
        return (1.0 / 6.0) * mV.Determinant();
    }

    /// <summary>
    /// Computes and returns the deformed node matrix.
    /// </summary>
    /// <returns>The deformed node matrix</returns>
    public MatrixXD GetNodeMatrix_Deformed()
    {
        MatrixXD Nx = new DenseMatrixXD(3, 4);
        Nx.SetColumn(0, Utils.ToVectorXD(Nodes[0].Pos));
        Nx.SetColumn(1, Utils.ToVectorXD(Nodes[1].Pos));
        Nx.SetColumn(2, Utils.ToVectorXD(Nodes[2].Pos));
        Nx.SetColumn(3, Utils.ToVectorXD(Nodes[3].Pos));
        return Nx;
    }

    /// <summary>
    /// Computes and returns the rest node matrix.
    /// </summary>
    /// <returns>The rest node matrix</returns>
    public MatrixXD GetNodeMatrix_Undeformed()
    {
        MatrixXD N0 = new DenseMatrixXD(3, 4);
        N0.SetColumn(0, Utils.ToVectorXD(Nodes[0].Pos0));
        N0.SetColumn(1, Utils.ToVectorXD(Nodes[1].Pos0));
        N0.SetColumn(2, Utils.ToVectorXD(Nodes[2].Pos0));
        N0.SetColumn(3, Utils.ToVectorXD(Nodes[3].Pos0));
        return N0;
    }

    /// <summary>
    /// Computes and returns the weights of the shape function S(p)
    /// at a point p inside the tetrahedron. The point p is expressed
    /// in iso-parametric coordinates.
    /// </summary>
    /// <param name="p">The point in iso-parametric coordinates</param>
    /// <returns>The interpolation weights of the shape function</returns>
    public VectorXD ComputeShapeFunctionValue(VectorXD p)
    {
        VectorXD S = new DenseVectorXD(4);
        S[0] = 1 - p[0] - p[1] - p[2];
        S[1] = p[0];
        S[2] = p[1];
        S[3] = p[2];
        return S;
    }

    /// <summary>
    /// Computes and returns the derivatives of the shape function
    /// w.r.t. a point in iso-parametric coordinates, G = dSdp.
    /// </summary>
    /// <returns>The matrix containing dSdp</returns>
    public MatrixXD ComputeShapeFunctionGradient()
    {
        MatrixXD G = new DenseMatrixXD(4, 3);
        G.Clear();
        G[0, 0] = -1;
        G[0, 1] = -1;
        G[0, 2] = -1;
        G[1, 0] = 1;
        G[2, 1] = 1;
        G[3, 2] = 1;
        return G;
    }

    /// <summary>
    /// Computes the part of the deformation gradient that
    /// only depends on the rest configuration and stores
    /// it in the the variable FRest.
    /// </summary>
    public void PrecomputeFRest()
    {
        MatrixXD N0 = this.GetNodeMatrix_Undeformed();
        MatrixXD G = this.ComputeShapeFunctionGradient();

        this.FRestToMatrix = G * (N0 * G).Inverse();
        this.FRestToVector = new DenseMatrixXD(9, 12);
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    FRestToVector[3 * j + k, 3 * i + k] = FRestToMatrix[i, j];
                }
            }
        }
    }

    #endregion

    #region MaterialAndDeformation

    /// <summary>
    /// Computes and returns the Cauchy strain E. 
    /// </summary>
    /// <returns>The Green strain tensor</returns>
    public MatrixXD ComputeCauchyStrain()
    {
        MatrixXD F = null;
        if (this.Simul.Corotational)
            F = ComputeCoRotDeformationGradient();
        else F = ComputeDefaultDeformationGradient();

        // What should I change here?

        return (0.5 * (F.Transpose() + F) - DenseMatrixXD.CreateIdentity(3));
    }

    /// <summary>
    /// Computes and returns the Green strain E. 
    /// </summary>
    /// <returns>The Green strain tensor</returns>
    public MatrixXD ComputeGreenStrain()
    {
        MatrixXD F = null;
        if (this.Simul.Corotational)
            F = ComputeCoRotDeformationGradient();
        else F = ComputeDefaultDeformationGradient();

        // What should I change here?

        return 0.5 * (F.Transpose() * F - DenseMatrixXD.CreateIdentity(3));
    }

    /// <summary>
    /// Computes and returns the deformation gradient F.
    /// </summary>
    /// <returns>The deformation gradient tensor</returns>
    public MatrixXD ComputeDefaultDeformationGradient()
    {
        MatrixXD Nx = GetNodeMatrix_Deformed();

        if (!Simul.Plastic) return Nx * this.FRestToMatrix;
        else return Nx * this.FRestToMatrix - FPlastic;
    }

    /// <summary>
    /// Computes and returns the deformation gradient F.
    /// </summary>
    /// <returns>The deformation gradient tensor</returns>
    public MatrixXD ComputeCoRotDeformationGradient()
    {
        MatrixXD Nx = GetNodeMatrix_Deformed();
        Nx = this.RotationOne.Transpose() * Nx;

        if (!Simul.Plastic) return Nx * this.FRestToMatrix;
        else return Nx * this.FRestToMatrix - FPlastic;
    }

    /// <summary>
    /// Computes and returns the first Lamé parameter.
    /// </summary>
    /// <returns>The first Lamé parameter</returns>
    public double ComputeLameParameterFirst()
    {
        return Simul.Young / (2 * (1 + Simul.Poisson));
    }

    /// <summary>
    /// Computes and returns the second Lamé parameter.
    /// </summary>
    /// <returns>The second Lamé parameter</returns>
    public double ComputeLameParameterSecond()
    {
        return (Simul.Young * Simul.Poisson) / ((1 + Simul.Poisson) * (1 - 2 * Simul.Poisson));
    }

    #endregion

    #region UpdateState

    /// <summary>
    /// Computes the rotation that best approximates the current state of
    /// the tetrahedron and stores it in the class variable RotationOne.
    /// Then, the diagonal 12x12 block matrix RotationAll is constructed
    /// to operate with the concatenation of tetrahedron nodes.
    /// </summary>
    public void UpdateRotation()
    {
        MatrixXD F = this.ComputeDefaultDeformationGradient();

        // Singular Value Decomposition: Decompose the deformation gradient in 
        // F = U*S*VT where U and and VT are orthogonal (rotation) matrix and S 
        // is a diagonal (scale) matrix. By recomputing F' = U*VT we are
        // elimitating the scale part.

        //Svd<double> factorSVD = F.Svd(true);
        //RotationOne = factorSVD.U * factorSVD.VT;

        // Gram-Schmidt orthogonalization: Computes the QR decomposition of the matrix
        // F = QR, where Q is an orthogonal (rotation) matrix and R is an upper triangular 
        // matrix. By keeping just Q we obtain an estimation of the rotation.

        GramSchmidt<double> factorGM = F.GramSchmidt();
        RotationOne = factorGM.Q;

        // Build the big block-diagonal rotation matrix.

        this.RotationAll = DenseMatrixXD.Create(12, 12, 0);
        for (int i = 0; i < 4; ++i)
            RotationAll.SetSubMatrix(3 * i, 3 * i, RotationOne);
    }

    /// <summary>
    /// Computes the current elastic strain. If the norm of the strain is higher
    /// than the YieldTol, a Creep part of the deformation must be accumulated
    /// as plastic deformation. If the norm of the plastic deformation exceeds 
    /// YieldMax, then the plastic deformation should be cap.
    /// </summary>
    public void UpdatePlasticity()
    {
        MatrixXD FElastic;
        if (Simul.Corotational)
            FElastic = this.ComputeCoRotDeformationGradient() - DenseMatrixXD.CreateIdentity(3);
        else FElastic = this.ComputeDefaultDeformationGradient() - DenseMatrixXD.CreateIdentity(3);

        // Conver part of the elastic deformation in plastic

        if (Simul.Plastic && FElastic.L2Norm() > Simul.YieldTol)
        {
            FPlastic = FPlastic + Simul.Creep * FElastic;
        }

        // Limit the maximum amount of plastic deformation

        if (Simul.Plastic && FPlastic.L2Norm() > Simul.YieldMax)
        {
            FPlastic = FPlastic * Simul.YieldMax / FPlastic.L2Norm();
        }
    }

    #endregion

    #region IsoparametricCoordinates

    /// <summary>
    /// Checks if the point passed as a parameter is
    /// inside the tetrahedron. The point is passed
    /// iso-parametric coordinates.
    /// </summary>
    /// <param name="p">The point in iso-parametric coordinates</param>
    /// <returns>True, if the point is inside the tetrahedron. False otherwise.</returns>
    public bool IsWithinElement(VectorXD p)
    {
        VectorXD S = this.ComputeShapeFunctionValue(p);

        float sum = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (S[i] < 0 - 1e-6 ||
                S[i] > 1 + 1e-6)
                return false;
            sum += (float)S[i];
        }

        return Mathf.Abs(sum - 1) < 1e-6;
    }

    /// <summary>
    /// Transforms the input point in rest coordinates to iso-parametric coordinates.
    /// </summary>
    /// <param name="x">The point to transform</param>
    /// <returns>The iso-parametric coordinates</returns>
    public VectorXD TransformRest2Isoparametric(VectorXD x)
    {
        return Rest2Iso * (x - Utils.ToVectorXD(Nodes[0].Pos0));
    }

    /// <summary>
    /// Transforms the input point in iso-parametric coordinates to rest coordinates.
    /// </summary>
    /// <param name="x">The point to transform</param>
    /// <returns>The rest coordinates</returns>
    public VectorXD TransformIsoparametric2Rest(VectorXD x)
    {
        return (Iso2Rest * x) + Utils.ToVectorXD(Nodes[0].Pos0);
    }

    /// <summary>
    /// Computes the change of basis matrix from rest coordinates to iso-parametric coordinates.
    /// </summary>
    /// <returns></returns>
    public MatrixXD ComputeRest2Isoparametric()
    {
        return this.ComputeIsoparametric2Rest().Inverse();
    }

    /// <summary>
    /// Computes the change of basis matrix from iso-parametric coordinates to rest coordinates.
    /// </summary>
    /// <returns></returns>
    public MatrixXD ComputeIsoparametric2Rest()
    {
        MatrixXD M = new DenseMatrixXD(3, 3);
        VectorXD dir01 = Utils.ToVectorXD(Nodes[1].Pos0 - Nodes[0].Pos0);
        VectorXD dir02 = Utils.ToVectorXD(Nodes[2].Pos0 - Nodes[0].Pos0);
        VectorXD dir03 = Utils.ToVectorXD(Nodes[3].Pos0 - Nodes[0].Pos0);
        M.SetColumn(0, dir01);
        M.SetColumn(1, dir02);
        M.SetColumn(2, dir03);
        return M;
    }

    #endregion

    #region LocalForcesAndJacobianAssembly

    /// <summary>
    /// Computes and returns the strain energy of the deformed tetrahedron.
    /// </summary>
    /// <returns>Return </returns>
    public double GetEnergy()
    {
        return this.ComputeEnergy();
    }

    /// <summary>
    /// Computes and assembles the elastic forces produced by this tetrahedron. 
    /// </summary>
    /// <param name="force">The globlal force vector where to assemble local force</param>
    public void GetForce(VectorXD force)
    {
        VectorXD tetForce = this.ComputeForce();

        for (int i = 0; i < 4; ++i)
            force.SetSubVector(Nodes[i].index, 3, force.SubVector(Nodes[i].index, 3) + tetForce.SubVector(3 * i, 3));
    }

    /// <summary>
    /// Computes and assembles the Jacobians of the elastic forces produced by this tetrahedron.
    /// </summary>
    /// <param name="DfDx">The global Jacobian matrix where to assemble local DfDx</param>
    /// <param name="DfDv">The global Jacobian matrix where to assemble local DfDv</param>
    public void GetJacobian(MatrixXD DfDx, MatrixXD DfDv)
    {
        MatrixXD tetDfDx = this.ComputeDfDx();
        MatrixXD tetDfDv = this.ComputeDfDv();

        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            {
                DfDx.SetSubMatrix(Nodes[i].index, Nodes[j].index, DfDx.SubMatrix(Nodes[i].index, 3, Nodes[j].index, 3) + tetDfDx.SubMatrix(3 * i, 3, 3 * j, 3));
                DfDv.SetSubMatrix(Nodes[i].index, Nodes[j].index, DfDv.SubMatrix(Nodes[i].index, 3, Nodes[j].index, 3) + tetDfDv.SubMatrix(3 * i, 3, 3 * j, 3));
            }
    }

    /// <summary>
    /// Computes and assembles the Jacobians of the elastic forces produced by this tetrahedron.
    /// </summary>
    /// <param name="DfDx">The global Jacobian matrix where to assemble local DfDx</param>
    /// <param name="DfDv">The global Jacobian matrix where to assemble local DfDv</param>
    public void GetJacobian(List<Triplet> DfDx, List<Triplet> DfDv)
    {
        MatrixXD tetDfDx = this.ComputeDfDx();
        MatrixXD tetDfDv = this.ComputeDfDv();

        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int ii = 0; ii < 3; ++ii)
                    for (int jj = 0; jj < 3; ++jj)
                    {
                        DfDx.Add(new Triplet(Nodes[i].index + ii, Nodes[j].index + jj, tetDfDx[3 * i + ii, 3 * j + jj]));
                        DfDv.Add(new Triplet(Nodes[i].index + ii, Nodes[j].index + jj, tetDfDv[3 * i + ii, 3 * j + jj]));
                    }
    }

    #endregion

    #region LocalEnergiesForcesAndJacobians

    /// <summary>
    /// Computes and returns the strain energy of the tetrahedron
    /// following the linear elastic model. NOTE: This refers to
    /// the total energy and not the energy density.
    /// </summary>
    /// <returns>The total strain energy of the tetrahedron</returns>
    public double ComputeEnergy()
    {

        MatrixXD e = this.ComputeCauchyStrain();

        return Volume0 * (firstlame * (e.Transpose() * e).Trace() + (secondlame / 2) * e.Trace() * e.Trace());
    }

    /// <summary>
    /// Computes and returns the forces at the nodes of the tetrahedron
    /// following the linear elastic model. All the forces concatenated
    /// in the same vector.
    /// </summary>
    /// <returns>The 12x1 vector of nodal forces</returns>
    public VectorXD ComputeForce()
    {

        // Compute stress density

        MatrixXD e = this.ComputeCauchyStrain();

        MatrixXD P = 2 * firstlame * e + secondlame * e.Trace() * DenseMatrixXD.CreateIdentity(3);

        // Compute force vector

        VectorXD p = Utils.ToVector_ColWise(P);

        VectorXD force = -Volume0 * p * FRestToVector;

        if (Simul.Corotational)
        {
            force = force * RotationAll.Transpose();
        }

        force += this.ComputeDamping();

        return force;
    }

    /// <summary>
    /// Computes and returns the damping force of this tetrahedron 
    /// considering it might be using the Rayleight damping model.
    /// </summary>
    /// <returns>The 12x1 damping force vector</returns>
    public VectorXD ComputeDamping()
    {
        if (Simul.AlphaDamping > 0 || Simul.BetaDamping > 0)
        {
            MatrixXD D = this.ComputeDfDv();
            VectorXD v = new DenseVectorXD(12);
            for (int i = 0; i < 0; ++i)
                v.SetSubVector(3 * i, 3, Utils.ToVectorXD(Nodes[i].Vel - Nodes[i].Pos0));

            return -D * v;
        }

        return DenseVectorXD.Create(12, 0);
    }

    /// <summary>
    /// Computes and returns the Jacobian of the forces at the nodes of
    /// the tetrahedron w.r.t. positions, following the linear elastic model. 
    /// All the forces are considered to be concatenated in the same vector.
    /// The derivative is computed w.r.t. the concatenated positions vector.
    /// </summary>
    /// <returns>The 12x12 matrix of force Jacobians w.r.t positions</returns>
    public MatrixXD ComputeDfDx()
    {

        MatrixXD dPdF = DenseMatrixXD.Create(9, 9, 0);
        dPdF[0, 0] = dPdF[4, 4] = dPdF[8, 8] = 2 * firstlame + secondlame;
        dPdF[4, 0] = dPdF[8, 0] = secondlame;
        dPdF[0, 4] = dPdF[8, 4] = secondlame;
        dPdF[0, 8] = dPdF[4, 8] = secondlame;

        MatrixXD muX = DenseMatrixXD.Create(3, 3, 0);
        muX[0, 0] = muX[0, 2] = firstlame;
        muX[2, 0] = muX[2, 2] = firstlame;
        muX[1, 1] = firstlame;
        dPdF.SetSubMatrix(1, 1, muX);
        dPdF.SetSubMatrix(5, 5, muX);
        dPdF[6, 2] = firstlame;
        dPdF[2, 6] = firstlame;

        MatrixXD jacobian = -Volume0 * this.FRestToVector.Transpose() * dPdF * this.FRestToVector;

        if (Simul.Corotational)
        {
            jacobian = RotationAll * jacobian * RotationAll.Transpose();
        }

        return jacobian;
    }

    /// <summary>
    /// Computes and returns the Jacobian of the forces at the nodes of
    /// the tetrahedron w.r.t. velocities, considering the FEM might be
    /// using the Rayleight damping model.
    /// </summary>
    /// <returns>The 12x12 matrix of force Jacobians w.r.t velocites</returns>
    public MatrixXD ComputeDfDv()
    {
        if (Simul.AlphaDamping > 0 || Simul.BetaDamping > 0)
        {
            MatrixXD J = this.ComputeDfDx();
            float mass = (float)this.Volume0 * Simul.Density / 4f;
            MatrixXD M = DenseMatrixXD.CreateDiagonal(12, 12, mass);
            MatrixXD D = Simul.AlphaDamping * M - Simul.BetaDamping * J;
            return D;
        }

        return DenseMatrixXD.Create(12, 12, 0);
    }

    #endregion

    #region Embedding

    /// <summary>
    /// If the tetrahedron has any embedded point from a mesh being used for visualization 
    /// or contact purposes (in embeddedPoints list) this method updates the position of 
    /// such point (embeddedPoint.position) through interpolation.
    /// </summary>
    public void UpdateEmbedded()
    {
        for (int i = 0; i < this.embeddedPoints.Count; ++i)
        {
            // COMPLETAR Parte 2

            VectorXD pos = new DenseVectorXD(3);
            VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);

            for(int j = 0; j < this.Nodes.Count; j++)
            {
                pos += Utils.ToVectorXD(Nodes[j].Pos) * shape[j];
            }
            embeddedPoints[i].position = pos-Utils.ToVectorXD(transform.root.position);
        }
    }

    /// <summary>
    /// Checks if the position emb.position is within the rest tetrahedron. The point
    /// is provided in rest coordinates. If the point is within the tetrahedron, sets 
    /// emb.pointIso to the iso-parametric coordinate of that point and stores the 
    /// embedded point in the list of embedded points this.embeddedPoints.
    /// </summary>
    /// <param name="emb">The embedded point</param>
    /// <returns>Returns true if the point is embedded</returns>
    public bool TryAddEmbedded(FEMSystem.EmbeddedPoint emb)
    {
        VectorXD position0 = emb.position;

        // COMPLETAR Parte 2
        VectorXD iso= Rest2Iso* (position0-Utils.ToVectorXD(Nodes[0].Pos0));
        if (IsWithinElement(iso))
        {
            emb.isoParam = iso;
            this.embeddedPoints.Add(emb);
            return true;
        }

        return false;
    }

    /// <summary>
    /// Adds to the total vector of forces the ones corresponding to the collision of 
    /// embedded points. For some given collision plane defined by a point p and a 
    /// normal n, this method checks if any of the points embedded in this tetrahedron
    /// (stored in embeddedPoints), reacts to the collision. In such case, it
    /// adds the force corresponding to that collision to the global force 
    /// vector. 
    /// </summary>
    /// <param name="force">The global force vector</param>
    /// <param name="k">Collision stiffness constant</param>
    /// <param name="p">The point defining the plane</param>
    /// <param name="n">The normal defining the plane</param>
    public void AddCollisionForceEmbedded(VectorXD force, double k, VectorXD p, VectorXD n)
    {
        for (int i = 0; i < embeddedPoints.Count; ++i)
        {
            // COMPLETAR Parte 3
            if ((embeddedPoints[i].position-p) * n < 0)
            {
                double d = (p-embeddedPoints[i].position) * n;
                VectorXD fp = k * n * d;
                VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);
                for(int j = 0; j < Nodes.Count; j++)
                {
                    int index = Nodes[j].index;
                    force.SetSubVector(index, 3, force.SubVector(index, 3) +fp*shape[j]);
                }
            }
        }
    }

    public void AddCollisionForceEmbeddedSphere(VectorXD force, double k, VectorXD p, double r)
    {
        for (int i = 0; i < embeddedPoints.Count; ++i)
        {
            // COMPLETAR Parte 3
            double norm = (embeddedPoints[i].position - p).L2Norm();
            if (norm < r)
            {
                VectorXD diff = embeddedPoints[i].position - p;
                VectorXD n = diff / norm;
                VectorXD fp = k * n * (r - norm);
                VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);
                for (int j = 0; j < Nodes.Count; j++)
                {
                    int index = Nodes[j].index;
                    force.SetSubVector(index, 3, force.SubVector(index, 3) + fp * shape[j]);
                }
            }
        }
    }

    /// <summary>
    /// Adds to the total force Jacobians the contribution corresponding to the collision of 
    /// embedded points. For some given collision plane defined by a point p and a normal n,
    /// this method checks if any of the points embedded in this tetrahedron (stored in 
    /// embeddedPoints), reacts to the collision. In such case, it adds the corresponding 
    /// Jacobian of the collision force to the global Jacobian 
    /// vector. 
    /// </summary>
    /// <param name="dfdx">The global dfdx matrix</param>
    /// <param name="dfdv">The global dfdv matrix</param>
    /// <param name="k">Collision stiffness constant</param>
    /// <param name="p">The point defining the plane</param>
    /// <param name="n">The normal defining the plane</param>
    public void AddCollisionJacobianEmbedded(MatrixXD dfdx, MatrixXD dfdv, double k, VectorXD p, VectorXD n)
    {
        MatrixXD block = -k * n.OuterProduct(n);
        for (int i = 0; i < embeddedPoints.Count; ++i)
        {
            // COMPLETAR Parte 3
            VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);
            if ((embeddedPoints[i].position - p) * n < 0)
            {
                for (int j = 0; j < Nodes.Count; j++)
                {
                    for (int ii = 0; ii < Nodes.Count; ii++)
                    {
                        MatrixXD blockij = block * shape[j] * shape[ii];

                        dfdx.SetSubMatrix(Nodes[j].index, Nodes[ii].index,
                            dfdx.SubMatrix(Nodes[j].index, 3, Nodes[ii].index, 3)
                            + blockij);

                        dfdx.SetSubMatrix(Nodes[ii].index, Nodes[j].index,
                            dfdx.SubMatrix(Nodes[ii].index, 3, Nodes[j].index, 3)
                            + blockij);
                    }
                }
            }
                
        }
    }

    public void AddCollisionJacobianEmbeddedSphere(MatrixXD dfdx, MatrixXD dfdv, double k, VectorXD p, double r)
    {
        
        for (int i = 0; i < embeddedPoints.Count; ++i)
        {
            // COMPLETAR Parte 3
            VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);
            VectorXD diff = embeddedPoints[i].position - p;
            double norm = diff.L2Norm();
            
            if (norm < r)
            {
                VectorXD n = diff / norm;
                MatrixXD block = -k * n.OuterProduct(n);
                for (int j = 0; j < Nodes.Count; j++)
                {
                    for (int ii = 0; ii < Nodes.Count; ii++)
                    {
                        MatrixXD blockij = block * shape[j] * shape[ii];
                        
                        dfdx.SetSubMatrix(Nodes[j].index, Nodes[ii].index,
                            dfdx.SubMatrix(Nodes[j].index, 3, Nodes[ii].index, 3)
                            + blockij);

                        dfdx.SetSubMatrix(Nodes[ii].index, Nodes[j].index,
                            dfdx.SubMatrix(Nodes[ii].index, 3, Nodes[j].index, 3)
                            + blockij);
                    }
                }
            }

        }
    }

    /// <summary>
    /// Adds to the total force Jacobians the contribution corresponding to the collision of 
    /// embedded points. For some given collision plane defined by a point p and a normal n,
    /// this method checks if any of the points embedded in this tetrahedron (stored in 
    /// embeddedPoints), reacts to the collision. In such case, it adds the corresponding 
    /// Jacobian of the collision force to the global Jacobian 
    /// vector. This is the sparse matrix version.
    /// </summary>
    /// <param name="dfdx">The global dfdx matrix (list of triplets)</param>
    /// <param name="dfdv">The global dfdv matrix (list of triplets)</param>
    /// <param name="k">Collision stiffness constant</param>
    /// <param name="p">The point defining the plane</param>
    /// <param name="n">The normal defining the plane</param>
    public void AddCollisionJacobianEmbedded(List<Triplet> dfdx, List<Triplet> dfdv, double k, VectorXD p, VectorXD n)
    {
        MatrixXD block = -k * n.OuterProduct(n);
        for (int i = 0; i < embeddedPoints.Count; ++i)
        {
            // OPCIONAL Parte 3
            VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);
            if ((embeddedPoints[i].position - p) * n < 0)
            {
                for (int j = 0; j < Nodes.Count; j++)
                {
                    for (int ii = 0; ii < Nodes.Count; ii++)
                    {
                        MatrixXD blockij = block * shape[j] * shape[ii];
                        for (int kk = 0; kk < 3; kk++)
                        {
                            for (int jj = 0; jj < 3; jj++)
                            {
                                dfdx.Add(new Triplet(Nodes[j].index+kk, Nodes[ii].index+jj, blockij[kk, jj]));
                                dfdx.Add(new Triplet(Nodes[ii].index+kk, Nodes[j].index+jj, blockij[kk, jj]));
                            }
                        }
                    }
                }
            }

        }
    }

    public void AddCollisionJacobianEmbeddedSphere(List<Triplet> dfdx, List<Triplet> dfdv, double k, VectorXD p, double r)
    {
        for (int i = 0; i < embeddedPoints.Count; ++i)
        {
            // OPCIONAL Parte 3
            VectorXD shape = ComputeShapeFunctionValue(embeddedPoints[i].isoParam);
            double norm = (embeddedPoints[i].position - p).L2Norm();
            VectorXD diff = embeddedPoints[i].position - p;
            if (norm < r)
            {
                VectorXD n = diff / norm;
                MatrixXD block = -k * n.OuterProduct(n);
                for (int j = 0; j < Nodes.Count; j++)
                {
                    for (int ii = 0; ii < Nodes.Count; ii++)
                    {
                        MatrixXD blockij = block * shape[j] * shape[ii];

                        for (int kk = 0; kk < 3; kk++)
                        {
                            for (int jj = 0; jj < 3; jj++)
                            {
                                dfdx.Add(new Triplet(Nodes[j].index + kk, Nodes[ii].index + jj, blockij[kk, jj]));
                                dfdx.Add(new Triplet(Nodes[ii].index + kk, Nodes[j].index + jj, blockij[kk, jj]));
                            }
                        }
                    }
                }
            }

        }
    }

    #endregion

    #region Testing

    /// <summary>
    /// Test the implementation of the current
    /// forces using finite-differences. Use this
    /// for checking correctness.
    /// </summary>
    public void TestForce()
    {
        VectorXD tetForce = this.ComputeForce();
        VectorXD tetForceFD = this.ComputeForce_FD();

        double normFD = tetForceFD.L2Norm();
        if (normFD < 1e-9)
            Debug.Log("Force error: " + normFD);
        else
        {
            VectorXD tetForceD = tetForce - tetForceFD;
            double error = tetForceD.L2Norm() / normFD;
            Debug.Log("Force error: " + error);
        }
    }

    /// <summary>
    /// Test the implementation of the current
    /// DfDx using finite-differences. Use this
    /// for checking correctness.
    /// </summary>
    public void TestDfDx()
    {
        MatrixXD tetDfDx = this.ComputeDfDx();
        MatrixXD tetDfDxFD = this.ComputeDfDx_FD();

        double normFD = tetDfDxFD.L2Norm();
        if (normFD < 1e-9)
            Debug.Log("DfDx error: " + normFD);
        else
        {
            MatrixXD tetDfDxD = tetDfDx - tetDfDxFD;
            double error = tetDfDxD.L2Norm() / normFD;
            Debug.Log("DfDx error: " + error);
        }
    }

    /// <summary>
    /// Computes the current force of the tetrahedron using finite-differences.
    /// </summary>
    /// <returns></returns>
    public VectorXD ComputeForce_FD()
    {
        float eps = 0.001f;

        VectorXD forcesFD = new DenseVectorXD(12);

        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                Nodes[i].Pos[j] += eps;

                double EnergyP = this.ComputeEnergy();

                Nodes[i].Pos[j] -= 2 * eps;

                double EnergyM = this.ComputeEnergy();

                Nodes[i].Pos[j] += eps;

                forcesFD[3 * i + j] = (EnergyM - EnergyP) / (2 * eps);
            }
        }

        return forcesFD;
    }

    /// <summary>
    /// Computes the current DfDx of the tetrahedron using finite-differences.
    /// </summary>
    /// <returns></returns>
    public MatrixXD ComputeDfDx_FD()
    {
        float eps = 0.001f;

        MatrixXD DfDxFD = new DenseMatrixXD(12);

        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                Nodes[i].Pos[j] += eps;

                VectorXD EnergyP = this.ComputeForce();

                Nodes[i].Pos[j] -= 2 * eps;

                VectorXD EnergyM = this.ComputeForce();

                Nodes[i].Pos[j] += eps;

                DfDxFD.SetColumn(3 * i + j, (EnergyP - EnergyM) / (2 * eps));
            }
        }

        return DfDxFD;
    }

    #endregion

    #region Other

    /// <summary>
    /// Paints the edges of the tetrahedron as Gizmos.
    /// </summary>
    private void OnDrawGizmos()
    {
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
            {
                Gizmos.color = Color.black;
                Gizmos.DrawLine(Nodes[i].Pos, Nodes[j].Pos);
            }
    }

    #endregion

}

