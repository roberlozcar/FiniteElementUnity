using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Fixer : MonoBehaviour
{

    // Use this for initialization
    void Start()
    {

    }

    // Update is called once per frame
    void Update()
    {

    }

    public bool IsInside(Vector3 pos)
    {
        return GetComponent<Collider>().bounds.Contains(pos);
    }
}
