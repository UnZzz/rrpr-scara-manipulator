using Godot;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System;

public partial class Arm : Node3D
{
    [ExportCategory("Joints")]
    [Export]
    public Node3D Joint1 { get; set; }
    [Export]
    public Node3D Joint2 { get; set; }
    [Export]
    public Node3D Joint3 { get; set; }
    [Export]
    public Node3D Joint4 { get; set; }
    [Export]
    public Node3D EndEffector { get; set; }

    [ExportCategory("Qs")]
    [Export]
    public float Theta1 { get; set; }
    [Export]
    public float Theta2 { get; set; }
    [Export]
    public float D3 { get; set; }
    [Export]
    public float Theta4 { get; set; }

    [ExportCategory("Parameters")]
    [Export]
    public float H { get; set; } = 8;
    [Export]
    public float A1 { get; set; } = 7;
    [Export]
    public float A2 { get; set; } = 14;
    [Export]
    public float D4 { get; set; } = 3;

    public DenseMatrix GetZRotateMatrix(float theta)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {Math.Cos(theta), -Math.Sin(theta), 0, 0},
            {Math.Sin(theta), Math.Cos(theta), 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetXRotateMatrix(float alpha)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, 0},
            {0, Math.Cos(alpha), -Math.Sin(alpha), 0},
            {0, Math.Sin(alpha), Math.Cos(alpha), 0},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetDTranslateMatrix(float d)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, d},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetATranslateMatrix(float a)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, a},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetHMatrix(float theta, float d, float a, float alpha)
    {
        var zRotate = GetZRotateMatrix(theta);
        var dTranslate = GetDTranslateMatrix(d);
        var aTranslate = GetATranslateMatrix(a);
        var xRotate = GetXRotateMatrix(alpha);

        var hMatrix = zRotate * dTranslate * aTranslate * xRotate;
        return hMatrix;
    }

    public void UpdateJointTransform(Node3D joint, DenseMatrix transform)
    {
        var basis = new Basis();
        
        var xAxis = new Vector3((float)transform[0, 0], (float)transform[1, 0], (float)transform[2, 0]);
        var yAxis = new Vector3((float)transform[0, 1], (float)transform[1, 1], (float)transform[2, 1]);
        var zAxis = new Vector3((float)transform[0, 2], (float)transform[1, 2], (float)transform[2, 2]);

        basis.X = xAxis; 
        basis.Y = yAxis;
        basis.Z = zAxis;

        var origin = new Vector3((float)transform[0, 3], (float)transform[1, 3], (float)transform[2, 3]);
        joint.GlobalTransform = new Transform3D(basis, origin);
    }

    public override void _Process(double delta)
    {
        var baseTransform = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        });

        var H01 = GetHMatrix(Theta1, 0, A1, 0);
        var joint1Transform = baseTransform * H01; 
        UpdateJointTransform(Joint1, joint1Transform);
        // GD.Print($"Joint 1 Transform: {joint1Transform}");

        var H12 = GetHMatrix(Theta2, 0, A2, 0);
        var joint2Transform = joint1Transform * H12;
        UpdateJointTransform(Joint2, joint2Transform);

        var H23 = GetHMatrix(0, -D3, 0, 0);
        var joint3Transform = joint2Transform * H23;
        UpdateJointTransform(Joint3, joint3Transform);

        var H34 = GetHMatrix(Theta4, -D4, 0, 0);
        var joint4Transform = joint3Transform * H34;
        UpdateJointTransform(Joint4, joint4Transform);

        // var H12 = GetHMatrix(Theta2, 0, A2, 0);
        // var joint2Transform = joint1Transform * H12;
        // UpdateJointTransform(Joint2, joint2Transform);

        // var H23 = GetHMatrix(0, -D3, 0, 0);
        // var joint3Transform = joint2Transform * H23;
        // UpdateJointTransform(Joint3, joint3Transform);

        // var H34 = GetHMatrix(Theta4, -D4, 0, 0);
        // var joint4Transform = joint3Transform * H34;
        // UpdateJointTransform(Joint4, joint4Transform);
    }
}
