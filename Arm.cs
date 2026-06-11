using System;
using Godot;
using MathNet.Numerics.Differentiation;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

[GlobalClass]
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
    public void SetTheta1(float value)
    {
        Theta1 = value;
    }
    [Export]
    public float Theta2 { get; set; }
    public void SetTheta2(float value)
    {
        Theta2 = value;
    }
    [Export]
    public float D3 { get; set; }
    public void SetD3(float value)
    {
        D3 = value;
    }
    [Export]
    public float Theta4 { get; set; }
    public void SetTheta4(float value)
    {
        Theta4 = value;
    }

    [ExportCategory("Parameters")]
    [Export]
    public float H { get; set; } = 8;
    [Export]
    public float A1 { get; set; } = 7;
    [Export]
    public float A2 { get; set; } = 14;
    [Export]
    public float D4 { get; set; } = 3;

    public DenseMatrix GetZRotateMatrix(double theta)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {Math.Cos(theta), -Math.Sin(theta), 0, 0},
            {Math.Sin(theta), Math.Cos(theta), 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetXRotateMatrix(double alpha)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, 0},
            {0, Math.Cos(alpha), -Math.Sin(alpha), 0},
            {0, Math.Sin(alpha), Math.Cos(alpha), 0},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetDTranslateMatrix(double d)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, d},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetATranslateMatrix(double a)
    {
        var mat = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, a},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        });
        return mat;
    }
    public DenseMatrix GetHMatrix(double theta, double d, double a, double alpha)
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
        var basis = new Basis
        {
            X = new Vector3((float)transform[0, 0], (float)transform[1, 0], (float)transform[2, 0]),
            Y = new Vector3((float)transform[0, 1], (float)transform[1, 1], (float)transform[2, 1]),
            Z = new Vector3((float)transform[0, 2], (float)transform[1, 2], (float)transform[2, 2])
        };

        var origin = new Vector3((float)transform[1, 3], (float)transform[2, 3], (float)transform[0, 3]);
        joint.GlobalTransform = new Transform3D(basis, origin);

        joint.RotateX(-Mathf.DegToRad(90));
        joint.RotateY(-Mathf.DegToRad(90));
    }

    // q=[theta1, theta2, d3, theta4]
    public DenseMatrix GetEETaskX(double[] q)
    {
        var h01 = GetHMatrix(0, H, 0, 0);
        var h12 = GetHMatrix(q[0], 0, A1, 0);
        var h23 = GetHMatrix(q[1], 0, A2, 0);
        var h34 = GetHMatrix(0, -q[2], 0, 0);
        var h4E = GetHMatrix(q[3], -D4, 0, 0);

        var EETransform = h01 * h12 * h23 * h34 * h4E;

        var xe = DenseMatrix.OfArray(new double[4, 1] {
            {EETransform[0, 3]}, //X
            {EETransform[1, 3]}, //Y
            {EETransform[2, 3]}, //Z
            {q[0] + q[1] + q[3]} //Phi
        });

        return xe;
    }

    private DenseMatrix GetJacobian(double[] q)
    {
        var Jacobian = new NumericalJacobian();

        Func<double[], double>[] f =
        {
            qq => GetEETaskX(qq)[0, 0], // x(q)
            qq => GetEETaskX(qq)[1, 0], // y(q)
            qq => GetEETaskX(qq)[2, 0], // z(q)
            qq => GetEETaskX(qq)[3, 0]  // phi(q)
        };

        var j = Jacobian.Evaluate(f, q);
        return DenseMatrix.OfArray(j);
    }

    public DenseMatrix GetCurrentJacobian()
    {
        var q = new double[] { Theta1, Theta2, D3, Theta4 };
        return GetJacobian(q);
    }

    public DenseMatrix GetCurrentEEX()
    {
        var q = new double[] { Theta1, Theta2, D3, Theta4 };
        return GetEETaskX(q);
    }

    public Vector3 GetCurrentEEPosition()
    {

        return EndEffector.GlobalTransform.Origin;
    }

    public override void _Process(double delta)
    {
        var baseTransform = DenseMatrix.OfArray(new double[4, 4] {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        });

        var H01 = GetHMatrix(0, H, 0, 0);
        var joint1Transform = baseTransform * H01;

        UpdateJointTransform(Joint1, joint1Transform);

        var H12 = GetHMatrix(Theta1, 0, A1, 0);
        var joint2Transform = joint1Transform * H12;
        UpdateJointTransform(Joint2, joint2Transform);

        var H23 = GetHMatrix(Theta2, 0, A2, 0);
        var joint3Transform = joint2Transform * H23;
        UpdateJointTransform(Joint3, joint3Transform);

        var H34 = GetHMatrix(0, -D3, 0, 0);
        var joint4Transform = joint3Transform * H34;
        UpdateJointTransform(Joint4, joint4Transform);

        var H4E = GetHMatrix(Theta4, -D4, 0, 0);
        var endEffectorTransform = joint4Transform * H4E;
        UpdateJointTransform(EndEffector, endEffectorTransform);
    }
}
