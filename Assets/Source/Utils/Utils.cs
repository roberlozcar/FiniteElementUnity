using UnityEngine;
using System.Collections;
using System.Collections.Generic;

using VectorXD = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MatrixXD = MathNet.Numerics.LinearAlgebra.Matrix<double>;
using DenseVectorXD = MathNet.Numerics.LinearAlgebra.Double.DenseVector;
using DenseMatrixXD = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix;
using Triplet = MathNet.Numerics.Tuple<int, int, double>;

static class Utils
{
    public static Vector3 ToVector3(VectorXD v)
    {
        Vector3 vOut = new Vector3();
        vOut.x = (float)v[0];
        vOut.y = (float)v[1];
        vOut.z = (float)v[2];
        return vOut;
    }

    public static VectorXD ToVectorXD(Vector3 v)
    {
        VectorXD vOut = new DenseVectorXD(3);
        vOut[0] = v.x;
        vOut[1] = v.y;
        vOut[2] = v.z;
        return vOut;
    }

    public static Vector3 computeHeatColorRGB(float alpha)
    {
        Vector3 heatHSL = computeHeatColorHSL(alpha);
        return convertHSLtoRGB(heatHSL); // Convert
    }

    public static Vector3 computeHeatColorHSL(float alpha)
    {
        float MAX = (240.0f / 360.0f);
        float h = (1.0f - alpha) * MAX;
        float s = 0.85f;
        float l = 0.5f;

        Vector3 hsl;
        hsl.x = h;
        hsl.y = s;
        hsl.z = l;
        return hsl;
    }

    public static Vector3 convertHSLtoRGB(Vector3 hsl)
    {
        double h = hsl[0];
        double sl = hsl[1];
        double l = hsl[2];

        double v;
        double r, g, b;

        r = l;   // default to gray
        g = l;
        b = l;
        v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
        if (v > 0)
        {
            double m;
            double sv;
            int sextant;
            double fract, vsf, mid1, mid2;

            m = l + l - v;
            sv = (v - m) / v;
            h *= 6.0;
            sextant = (int)h;
            fract = h - sextant;
            vsf = v * sv * fract;
            mid1 = m + vsf;
            mid2 = v - vsf;
            switch (sextant)
            {
                case 0:
                    r = v;
                    g = mid1;
                    b = m;
                    break;
                case 1:
                    r = mid2;
                    g = v;
                    b = m;
                    break;
                case 2:
                    r = m;
                    g = v;
                    b = mid1;
                    break;
                case 3:
                    r = m;
                    g = mid2;
                    b = v;
                    break;
                case 4:
                    r = mid1;
                    g = m;
                    b = v;
                    break;
                case 5:
                    r = v;
                    g = m;
                    b = mid2;
                    break;
            }
        }

        Vector3 rgb;
        rgb.x = (float)r;
        rgb.y = (float)g;
        rgb.z = (float)b;
        return rgb;
    }

    /// <summary>
    /// Transforms the input matrix to a vector column-wise.
    /// </summary>
    /// <param name="M">The matrix to transform</param>
    /// <returns>The column-wise matrix vector</returns>
    public static VectorXD ToVector_ColWise(MatrixXD M)
    {
        return new DenseVectorXD(M.ToColumnWiseArray());
    }

    /// <summary>
    /// Transforms the input matrix to a vector row-wise.
    /// </summary>
    /// <param name="M">The matrix to transform</param>
    /// <returns>The row-wise matrix vector</returns>
    public static VectorXD ToVector_RowWise(MatrixXD M)
    {
        return new DenseVectorXD(M.ToRowWiseArray());
    }


    public static List<Triplet> CollapseTriplets(int numRow, int numCol, List<Triplet> triplets)
    {
        Dictionary<int, double>[] rows = new Dictionary<int, double>[numRow];
        for (int i = 0; i < numRow; ++i)
            rows[i] = new Dictionary<int, double>();

        for (int i = 0; i < triplets.Count; ++i)
        {
            if (rows[triplets[i].Item1].ContainsKey(triplets[i].Item2))
                rows[triplets[i].Item1][triplets[i].Item2] += triplets[i].Item3;
            else rows[triplets[i].Item1][triplets[i].Item2] = triplets[i].Item3;
        }

        List<Triplet> collapsed = new List<Triplet>(triplets.Count);
        for (int i = 0; i < numRow; ++i)
            foreach (KeyValuePair<int, double> entry in rows[i])
            {
                collapsed.Add(new Triplet(i, entry.Key, entry.Value));
            }

        return collapsed;
    }

}

