using OBJLoaderPackage;
using System.IO;
using UnityEngine;

public class ObjFromFile : MonoBehaviour
{
    void Start()
    {
        //file path
        string filePath = @"I:\random\cylinder.obj";
        if (!File.Exists(filePath))
        {
            Debug.LogError("Please set FilePath in ObjFromFile.cs to a valid path.");
            return;
        }
        //create stream and load
        var loadedObj = new OBJLoader().Load(filePath);
    }
}
