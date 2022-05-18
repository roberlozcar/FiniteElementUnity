using System.Collections;
using System.Globalization;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System;

using OBJLoaderPackage;

public class ObjectLoader : MonoBehaviour
{
    public string Ident = "";

    [Header("Parameters")]
    public bool Corotational = false;  // Use co-rotational method
    public bool Plastic = false; // Use plasticity
    public float YieldTol = 0.01f; // Plasticity yield tolerance parameter
    public float YieldMax = 0.4f; // Plasticity yield maximum parameter
    public float Creep = 0.8f; // Plasticity creep parameters
    public float Density = 1000.0f; // Density value
    public float Young = 2.5e5f; // Young modulus
    public float Poisson = 0.35f; // Poisson's ratio
    public float Drag = 0.0f; // Air drag
    public float AlphaDamping = 0.2f; // Alpha parameter for the Rayleight damping
    public float BetaDamping = 0.001f; // Beta parameter for the Rayleight damping
    public FEMSystem.CollisionType CollType = FEMSystem.CollisionType.SimulatedNode; // Collision type
    public float CollStiffness = 1000000; // The collision forces stiffness
    public Vector3 Gravity = new Vector3(0.0f, -9.81f, 0.0f);

    [HideInInspector] public List<GameObject> CollObstacles; // Collection of obstacles (planes)
    [HideInInspector] public GameObject go;
    [HideInInspector] public List<Vector3> vSimNodes;
    [HideInInspector] public List<int[]> vSimIndices;
    [HideInInspector] public FEMSystem simulable = null;

    void Awake()
    {
        if (!this.ReadSimulationMeshes(Ident, out go, out vSimNodes, out vSimIndices))
        {
            Debug.LogError("Error reading meshes!: " + name);
        }

        go.transform.parent = this.transform;

        // Create the simulable and embedded the 
        // visualization mesh in the simulation

        Mesh mesh = go.transform.GetChild(0).GetComponent<MeshFilter>().mesh;
        simulable = CreateSimulable();
        simulable.SetEmbeddedMesh(mesh);
    }

    /// <summary>
    /// Given some resource name, this method loads from the Resources folder
    /// the visualization mesh ([name]_viz.obj) and the simulation mesh (stored in
    /// two files ([name]_sim.1.node and [name]_sim.1.ele).
    /// </summary>
    bool ReadSimulationMeshes(string name,
                              out GameObject go, 
                              out List<Vector3> vSimNodes, 
                              out List<int[]> vSimIndices)
    {
        // Initialize output parameters

        go = null;
        vSimNodes = new List<Vector3>();
        vSimIndices = new List<int[]>();

        ///////////////////////////////////////////////////
        // Read the visualization mesh and create object //
        ///////////////////////////////////////////////////

        string pathViz = "./Assets/Resources/" + name + "_viz.obj";
        if (!File.Exists(pathViz))
        {
            Debug.LogError("Resource not found!");
            return false;
        }

        // Load object from the .obj file

        go = new OBJLoader().Load(pathViz);
        go.transform.position += this.transform.position;
        ///////////////////////////////////////////////
        // Read the simulation mesh and create lists //
        ///////////////////////////////////////////////

        string pathSimNodes = "./Assets/Resources/" + name + "_sim.1.node";
        string pathSimTetras = "./Assets/Resources/" + name + "_sim.1.ele";
        if (!File.Exists(pathSimNodes))
        {
            Debug.LogError("Resource not found!: " + pathSimNodes);
            return false;
        }
        if (!File.Exists(pathSimTetras))
        {
            Debug.LogError("Resource not found!: " + pathSimTetras);
            return false;
        }

        string[] nodesContentLines = File.ReadAllLines(pathSimNodes);
        string[] tetrasContentLines = File.ReadAllLines(pathSimTetras);

        if (!ParseTetGenNodes(nodesContentLines, out vSimNodes))
        {
            Debug.LogError("Error reading nodes!: " + pathSimNodes);
            return false;
        }

        if (!ParseTetGenTetras(tetrasContentLines, out vSimIndices))
        {
            Debug.LogError("Error reading nodes!: " + pathSimNodes);
            return false;
        }

        return true;
    }

    #region CreationAndParsing

    /// <summary>
    /// Creates the simulable object. Takes the positions of the tetrahedral
    /// mesh vertices (vSimVertices) and the indices (vSimIndices) and creates
    /// the corresponding Node/Tetrahedron objects and components. 
    /// </summary>
    /// <returns>The created simulable</returns>
    FEMSystem CreateSimulable()
    {
        FEMSystem femSystem = this.gameObject.AddComponent<FEMSystem>();

        // Set parameters from editor

        // COMPLETAR Parte 1

        femSystem.Corotational = Corotational;
        femSystem.Plastic = Plastic;
        femSystem.YieldTol = YieldTol;
        femSystem.YieldMax = YieldMax;
        femSystem.Creep=Creep;
        femSystem.Density=Density;
        femSystem.Young=Young;
        femSystem.Poisson=Poisson;
        femSystem.Drag = Drag;
        femSystem.AlphaDamping = AlphaDamping;
        femSystem.BetaDamping = BetaDamping;
        femSystem.BetaDamping = BetaDamping;
        femSystem.CollStiffness = CollStiffness;
        femSystem.CollStiffness = CollStiffness;
        femSystem.Gravity = Gravity;
        femSystem.CollObstacles = CollObstacles;
        femSystem.CollType = CollType;


    // Create nodes from vSimNodes
        femSystem.Nodes = new List<Node>();

        // COMPLETAR Parte 1

        for (int i = 0; i < vSimNodes.Count; i++)
        {
            GameObject nodego = new GameObject();
            nodego.transform.parent = femSystem.transform;
            nodego.AddComponent<Node>();
            Node node = nodego.GetComponent<Node>();
            node.transform.localPosition = vSimNodes[i];
            node.Pos = vSimNodes[i];
            femSystem.Nodes.Add(node);
        }

        // Create tetrahedra from vSimIndices
		femSystem.Tetras = new List<Tetrahedron>();

        // COMPLETAR Parte 1
        
        for (int i = 0; i < vSimIndices.Count; i++)
        {
            GameObject tetrago = new GameObject();
            tetrago.AddComponent<Tetrahedron>();
            tetrago.transform.parent = femSystem.transform;
            Tetrahedron tetra = tetrago.GetComponent<Tetrahedron>();
            tetra.Nodes = new List<Node>();
            tetra.Nodes.Add(femSystem.Nodes[vSimIndices[i][0]]);
            tetra.Nodes.Add(femSystem.Nodes[vSimIndices[i][1]]);
            tetra.Nodes.Add(femSystem.Nodes[vSimIndices[i][2]]);
            tetra.Nodes.Add(femSystem.Nodes[vSimIndices[i][3]]);
            femSystem.Tetras.Add(tetra);
        }

        return femSystem;
    }

    /// <summary>
    /// This method parses the provided array of file lines to extract the 
    /// vertices of a tet-mesh. It returns true if everything went well.
    /// </summary>
    bool ParseTetGenNodes(string[] nodesContentLines, out List<Vector3> vNodes)
    {
        vNodes = new List<Vector3>();

		// NOTA: 
		// Para parsear números flotantes hay que tener en
		// cuenta el formato de número en el que está escrito.
		// Para ello hay que instanciar un objeto de la clase
		// CultureInfo que almacena información de localización. 
		// Los números con "." como separador decimal como 1.425
		// tienen localización de EEUU, "en-US".
		CultureInfo locale = new CultureInfo("en-US");

        // COMPLETAR Parte 1
        char[] separator = { ' ' };
        for(int i = 1; i < nodesContentLines.Length-1; i++)
        {
            String[] split = nodesContentLines[i].Split(separator, StringSplitOptions.RemoveEmptyEntries);
            float x = Single.Parse(split[1], locale);
            float y = Single.Parse(split[2], locale);
            float z = Single.Parse(split[3], locale);
            vNodes.Add(new Vector3(x,y,z));
        }

        return true;
    }

    /// <summary>
    /// This method parses the provided array of file lines to extract the 
    /// indices of a tet-mesh. It returns true if everything went well.
    /// NOTE: Consider TetGen stores nodes with base 1 not base 0.
    /// </summary>
    bool ParseTetGenTetras(string[] tetrasContentLines, out List<int[]> vTetras)
    {
        vTetras = new List<int[]>();

        // COMPLETAR Parte 1

        char[] separator = { ' ' };
        for (int i = 1; i<tetrasContentLines.Length - 1; i++)
        {
            String[] split = tetrasContentLines[i].Split(separator, 
                StringSplitOptions.RemoveEmptyEntries);
            vTetras.Add(new int[] {int.Parse(split[1])-1, int.Parse(split[2])-1, 
                int.Parse(split[3])-1, int.Parse(split[4])-1});
        }

        return true;
    }

    #endregion
}
