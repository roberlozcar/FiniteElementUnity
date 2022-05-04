using System.Collections.Generic;
using UnityEngine;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using SparseChol = MathNet.Numerics.LinearAlgebra.Double.Factorization.SparseChol;
using SparseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.SparseMatrix;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;
using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;

/// <summary>
/// Basic physics manager capable of simulating a given ISimulable
/// implementation using diverse integration methods: explicit,
/// implicit, Verlet and semi-implicit.
/// </summary>
public class PhysicsManager : MonoBehaviour 
{
    public enum SolverType
    {
        Dense,
        Sparse,
    }

	/// <summary>
	/// Default constructor. Zero all. 
	/// </summary>
	public PhysicsManager()
	{
		Paused = true;
		TimeStep = 0.01f;
		Gravity = new Vector3 (0.0f, -9.81f, 0.0f);
		IntegrationMethod = Integration.Symplectic;
	}

	/// <summary>
	/// Integration method.
	/// </summary>
	public enum Integration
	{
		Symplectic = 0,
        Implicit = 1,
	};

    #region InEditorVariables

    [Header("Objects")]
    public List<ObjectLoader> LoadedObjects;
    public List<GameObject> CollObstacles;
    public List<Fixer> NodeFixers;
    [Space(5)]

    [Header("Integration")]
    public bool Paused = true;
    public bool Step = false;
	public float TimeStep = 0.05f;
    public Integration IntegrationMethod = Integration.Implicit;
    public SolverType solverType = SolverType.Dense;
    [Space(5)]

    [Header("Debug grid")]
    public bool CreateGrid = false;
    public bool FixPoints = true;
    public Vector3Int Dims = new Vector3Int(3, 3, 3);
    public Vector3 Gravity = new Vector3(0, -9.81f, 0);
    public float Young = 2.5e5f;
    public float Poisson = 0.35f;
    public float Density = 1000.0f;
    public float Drag = 0.0f;
    public float AlfaDamping = 0.2f;
    public float BetaDamping = 0.001f;
    public bool Corotational = false;
    public bool Plastic = false;
    public float YieldTol = 0.01f;
    public float YieldMax = 0.4f;
    public float Creep = 0.8f;
    public float CollStifness = 10000000.0f;
    [Space(5)]

    #endregion

    #region OtherVariables

    private List<GameObject> m_simObjects = new List<GameObject>();
    private List<ISimulable> m_simulables = new List<ISimulable>();
    private int m_numDoFs;

    #endregion

    #region MonoBehaviour

    public void Start()
    {
        if (CreateGrid)
        {
            this.CreateDebugGrid();
        }
        else
        {
            for (int i = 0; i < this.LoadedObjects.Count; ++i)
                this.m_simObjects.Add(this.LoadedObjects[i].gameObject);
        }

        // Parse the simulable objects and initialize their state indices

        m_numDoFs = 0;

        foreach (GameObject obj in m_simObjects)
        {
            ISimulable simobj = obj.GetComponent<ISimulable>();
            if (simobj != null)
            {
                m_simulables.Add(simobj);

                // Initialize simulable object
                simobj.Initialize(m_numDoFs, this, NodeFixers, CollObstacles);

                // Retrieve pos and vel size
                m_numDoFs += simobj.GetNumDoFs();
            }
        }
    }

    public void Update()
    {
        if (Input.GetKeyUp(KeyCode.P))
            this.Paused = !this.Paused;

        if (Input.GetKeyUp(KeyCode.S))
            this.Step = !this.Step;

        if (this.Step)
            this.Paused = false;

        if (this.Paused)
            return; // Not simulating

        // Select integration method
        switch (this.IntegrationMethod)
        {
            case Integration.Symplectic: this.StepSymplectic(); break;
            case Integration.Implicit:
                switch (this.solverType)
                {
                    case SolverType.Dense: this.StepImplicitDense(); break;
                    case SolverType.Sparse: this.StepImplicitSparse(); break;
                }
                break;
            default:
                throw new System.Exception("[ERROR] Should never happen!");
        }

        if (this.Step)
        {
            this.Paused = true;
            this.Step = false;
        }
    }

    #endregion

    /// <summary>
    /// Performs a simulation step using Symplectic integration.
    /// </summary>
    private void StepSymplectic()
	{
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        f.Clear();
        MatrixXD Minv = new DenseMatrixXD(m_numDoFs);
        Minv.Clear();

        foreach (ISimulable obj in m_simulables)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetMassInverse(Minv);
        }

        foreach (ISimulable obj in m_simulables)
        {
            obj.FixVector(f);
        }

        v += TimeStep * (Minv * f);
        x += TimeStep * v;

        foreach (ISimulable obj in m_simulables)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Implicit integration and dense algebra matrices.
    /// </summary>
    private void StepImplicitDense()
    {
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        f.Clear();

        MatrixXD mass = new DenseMatrixXD(m_numDoFs);
        MatrixXD dfdx = new DenseMatrixXD(m_numDoFs);
        MatrixXD dfdv = new DenseMatrixXD(m_numDoFs);
        dfdx.Clear();
        dfdv.Clear();
        mass.Clear();

        foreach (ISimulable obj in m_simulables)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetJacobian(dfdx, dfdv);
            obj.GetMass(mass);
        }

        MatrixXD A = mass - TimeStep * dfdv - (TimeStep * TimeStep) * dfdx;
        VectorXD b = (mass - TimeStep * dfdv) * v + TimeStep * f;

        //MatrixXD test = new DenseMatrixXD(m_numDoFs);
        //A.CopyTo(test);
        //Debug.Log("A: " + test.L2Norm());
        //Debug.Log("b: " + b.L2Norm());

        foreach (ISimulable obj in m_simulables)
        {
            obj.FixMatrix(A);
            obj.FixVector(b);
        }

        v = A.Solve(b);
        x += TimeStep * v;

        foreach (ISimulable obj in m_simulables)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

    /// <summary>
    /// Performs a simulation step using Implicit integration and sparse
    /// algebra matrices using a basic a list of triplets and a Cholesky
    /// sparse solver provided by CSparse (theoretically the fastest
    /// option but nothing too crazy).
    /// </summary>
    private void StepImplicitSparse()
    {
        VectorXD x = new DenseVectorXD(m_numDoFs);
        VectorXD v = new DenseVectorXD(m_numDoFs);
        VectorXD f = new DenseVectorXD(m_numDoFs);
        f.Clear();
        List<Triplet> massTriplets = new List<Triplet>();
        List<Triplet> dfdxTriplets = new List<Triplet>();
        List<Triplet> dfdvTriplets = new List<Triplet>();

        bool[] stencil = new bool[m_numDoFs];
        for (int i = 0; i < m_numDoFs; ++i)
            stencil[i] = false;

        foreach (ISimulable obj in m_simulables)
        {
            obj.GetPosition(x);
            obj.GetVelocity(v);
            obj.GetForce(f);
            obj.GetFixedStencil(stencil);
            obj.GetJacobian(dfdxTriplets, dfdvTriplets);
            obj.GetMass(massTriplets);
        }

        SparseMatrixXD mass = this.BuildSparse(massTriplets, m_numDoFs, m_numDoFs);
        SparseMatrixXD dfdx = this.BuildSparse(dfdxTriplets, m_numDoFs, m_numDoFs);
        SparseMatrixXD dfdv = this.BuildSparse(dfdvTriplets, m_numDoFs, m_numDoFs);

        SparseMatrixXD A = mass - TimeStep * TimeStep * dfdx - TimeStep * dfdv;
        VectorXD b = (mass - TimeStep * dfdv) * v + TimeStep * f;

        //MatrixXD test = new DenseMatrixXD(m_numDoFs);
        //A.CopyTo(test);
        //Debug.Log("A: " + test.L2Norm());
        //Debug.Log("b: " + b.L2Norm());

        this.SetFixedInVector(stencil, b);
        this.SetFixedInMatrix(stencil, A);

        SparseChol solver = SparseChol.Create(A, CSparse.ColumnOrdering.Natural);

        v = solver.Solve(b);
        x += TimeStep * v;

        foreach (ISimulable obj in m_simulables)
        {
            obj.SetPosition(x);
            obj.SetVelocity(v);
        }
    }

    public SparseMatrixXD BuildSparse(List<Triplet> vc, int N, int M)
    {
        System.Collections.Generic.Dictionary<System.Tuple<int, int>, double> map
            = new System.Collections.Generic.Dictionary<System.Tuple<int, int>, double>();
        foreach (Triplet triplet in vc)
        {
            System.Tuple<int, int> key = new System.Tuple<int, int>(triplet.Item1, triplet.Item2);
            if (map.ContainsKey(key))
                map[key] += triplet.Item3;
            else map.Add(key, triplet.Item3);
        }

        vc.Clear();

        foreach (var element in map)
            vc.Add(new Triplet(element.Key.Item1, element.Key.Item2, element.Value));

        return SparseMatrixXD.OfIndexed(N, M, vc);
    }

    public void SetFixedInVector(bool[] stencil, VectorXD vector)
    {
        for (int i = 0; i < vector.Count; ++i)
            if (stencil[i])
                vector[i] = 0;
    }

    public void SetFixedInMatrix(bool[] stencil, MatrixXD matrix)
    {
        for (int i = 0; i < matrix.RowCount; ++i)
            for (int j = 0; j < matrix.ColumnCount; ++j)
                if (stencil[i] || stencil[j])
                    if (i == j)
                        matrix[i, j] = 1;
                    else matrix[i, j] = 0;
    }

    public void SetFixedInMatrix(bool[] stencil, SparseMatrixXD matrix)
    {
        matrix.MapIndexedInplace((i, j, x) => {
            if (stencil[i] || stencil[j])
            {
                if (i == j)
                    return 1;
                if (i != j)
                    return 0;
            }
            return x;
        }, MathNet.Numerics.LinearAlgebra.Zeros.AllowSkip);
    }

    void CreateDebugGrid()
    {
        GameObject femSystemGO = new GameObject();
        femSystemGO.AddComponent<FEMSystem>();
        femSystemGO.name = "DebugGrid";
        FEMSystem femSystem = femSystemGO.GetComponent<FEMSystem>();
        femSystem.Nodes = new List<Node>();
        femSystem.Tetras = new List<Tetrahedron>();
        femSystem.Young = Young;
        femSystem.Poisson = Poisson;
        femSystem.Density = Density;
        femSystem.Drag = Drag;
        femSystem.AlphaDamping = AlfaDamping;
        femSystem.BetaDamping = BetaDamping;
        femSystem.Corotational = Corotational;
        femSystem.Plastic = Plastic;
        femSystem.YieldMax = YieldMax;
        femSystem.YieldTol = YieldTol;
        femSystem.Creep = Creep;
        femSystem.Manager = this;
        femSystem.CollStiffness = CollStifness;

        femSystem.transform.parent = this.transform;

        this.m_simObjects = new List<GameObject>();
        this.m_simObjects.Add(femSystemGO);

        // Create nodes

        List<List<List<GameObject>>> mnodesGO = new List<List<List<GameObject>>>();

        for (int i = 0; i < Dims.x; ++i)
        {
            mnodesGO.Add(new List<List<GameObject>>());
            for (int j = 0; j < Dims.y; ++j)
            {
                mnodesGO[i].Add(new List<GameObject>());
                for (int k = 0; k < Dims.z; ++k)
                {
                    GameObject nodeGO = new GameObject();
                    nodeGO.transform.position = new Vector3(6*i, 6*j, 6*k);
                    nodeGO.name = "Node" + femSystem.Nodes.Count;
                    mnodesGO[i][j].Add(nodeGO);

                    nodeGO.AddComponent<Node>();
                    nodeGO.transform.parent = femSystemGO.transform;
                    femSystem.Nodes.Add(nodeGO.GetComponent<Node>());

                    if (i == 0 && FixPoints)
                        nodeGO.GetComponent<Node>().Fixed = true;
                }
            }
        }

        // Create tetrahedrons

        for (int k = 0; k < Dims.z - 1; ++k)
        {
            for (int j = 0; j < Dims.y - 1; ++j)
            {
                for (int i = 0; i < Dims.x - 1; ++i)
                {
                    List<Node> Nodes = new List<Node>();
                    Nodes.Add(mnodesGO[i][j][k].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i + 1][j][k].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i + 1][j + 1][k].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i][j + 1][k].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i][j][k + 1].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i + 1][j][k + 1].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i + 1][j + 1][k + 1].GetComponent<Node>());
                    Nodes.Add(mnodesGO[i][j + 1][k + 1].GetComponent<Node>());

                    List<List<int>> Tetras = new List<List<int>>(6);
                    for (int t = 0; t < 6; ++t)
                        Tetras.Add(new List<int>());
                    Tetras[0].Add(0); Tetras[0].Add(4); Tetras[0].Add(1); Tetras[0].Add(3);
                    Tetras[1].Add(3); Tetras[1].Add(4); Tetras[1].Add(1); Tetras[1].Add(2);
                    Tetras[2].Add(3); Tetras[2].Add(4); Tetras[2].Add(2); Tetras[2].Add(7);
                    Tetras[3].Add(7); Tetras[3].Add(4); Tetras[3].Add(2); Tetras[3].Add(6);
                    Tetras[4].Add(1); Tetras[4].Add(4); Tetras[4].Add(5); Tetras[4].Add(2);
                    Tetras[5].Add(5); Tetras[5].Add(2); Tetras[5].Add(4); Tetras[5].Add(6);

                    for (int t = 0; t < 6; ++t)
                    {
                        GameObject tetGO = new GameObject();
                        tetGO.AddComponent<Tetrahedron>();
                        tetGO.name = "Tetra" + femSystem.Tetras.Count;
                        tetGO.transform.parent = femSystemGO.transform;
                        Tetrahedron tet = tetGO.GetComponent<Tetrahedron>();
                        tet.Nodes = new List<Node>();
                        for (int n = 0; n < 4; ++n)
                            tet.Nodes.Add(Nodes[Tetras[t][n]]);

                        femSystem.Tetras.Add(tet);
                    }
                }
            }
        }
    }

}
