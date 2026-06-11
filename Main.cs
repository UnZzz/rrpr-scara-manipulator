using System;
using System.Reflection.PortableExecutable;
using Godot;
using MathNet.Numerics.Differentiation;
using MathNet.Numerics.LinearAlgebra.Double;

public enum IKMethod
{
	PseudoInverse,
	DampedLeastSquaresWithGain,
	JacobianTranspose
}
public partial class Main : Node
{

	[Export]
	public MultiMesh MultiMesh { get; set; }

	[Export]
	public IKMethod IKMethod { get; set; } = IKMethod.PseudoInverse;
	public void SetIKMethod(int value)
	{
		IKMethod = (IKMethod)value;
	}

	[Export]
	public double StepSize { get; set; } = 0.01;
	public void SetStepSize(double value)
	{
		StepSize = value;
	}


	[Export]
	public Arm Arm { get; set; }
	[Export]
	public Node3D Target { get; set; }

	[Export]
	public double TargetX { get; set; }
	public void SetTargetX(double value)
	{
		TargetX = value;
	}
	[Export]
	public double TargetY { get; set; }
	public void SetTargetY(double value)
	{
		TargetY = value;
	}
	[Export]
	public double TargetZ { get; set; }
	public void SetTargetZ(double value)
	{
		TargetZ = value;
	}
	[Export]
	public double TargetTheta { get; set; }
	public void SetTargetTheta(double value)
	{
		TargetTheta = value;
	}

	public void PseudoInverseIK(double stepSize = 1.0)
	{
		var jacobian = Arm.GetCurrentJacobian();
		var inverseJacobian = jacobian.PseudoInverse();

		var currentEEX = Arm.GetCurrentEEX();
		var targetEEX = DenseMatrix.OfArray(new double[4, 1] {
			{ TargetX },
			{ TargetY },
			{ TargetZ },
			{ TargetTheta }
		});

		var error = targetEEX - currentEEX;

		var deltaQ = inverseJacobian * error;

		Arm.SetTheta1(Arm.Theta1 + (float)(stepSize * deltaQ[0, 0]));
		Arm.SetTheta2(Arm.Theta2 + (float)(stepSize * deltaQ[1, 0]));
		Arm.SetD3(Arm.D3 + (float)(stepSize * deltaQ[2, 0]));
		Arm.SetTheta4(Arm.Theta4 + (float)(stepSize * deltaQ[3, 0]));
	}

	public void DampedLeastSquaresIKWithGain(double lambda = 0.1, double stepSize = 1.0)
	{
		var jacobian = Arm.GetCurrentJacobian();

		var currentEEX = Arm.GetCurrentEEX();

		var targetEEX = DenseMatrix.OfArray(new double[4, 1] {
		{ TargetX },
		{ TargetY },
		{ TargetZ },
		{ TargetTheta }
	});

		var error = targetEEX - currentEEX;

		var K = DenseMatrix.OfArray(new double[4, 4] {
		{ 1.0, 0.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 0.0, 1.0, 0.0 },
		{ 0.0, 0.0, 0.0, 1.0 }
	});

		var v = K * error;

		var JT = jacobian.Transpose();

		var identity = DenseMatrix.CreateIdentity(jacobian.RowCount);
		var A = jacobian * JT + lambda * lambda * identity;

		var temp = A.Solve(v);
		var deltaQ = JT * temp;

		Arm.SetTheta1(Arm.Theta1 + (float)(stepSize * deltaQ[0, 0]));
		Arm.SetTheta2(Arm.Theta2 + (float)(stepSize * deltaQ[1, 0]));
		Arm.SetD3(Arm.D3 + (float)(stepSize * deltaQ[2, 0]));
		Arm.SetTheta4(Arm.Theta4 + (float)(stepSize * deltaQ[3, 0]));
	}

	public void JacobianTransposeIK(double stepSize = 0.01)
	{
		var jacobian = Arm.GetCurrentJacobian();

		var currentEEX = Arm.GetCurrentEEX();

		var targetEEX = DenseMatrix.OfArray(new double[4, 1] {
		{ TargetX },
		{ TargetY },
		{ TargetZ },
		{ TargetTheta }
	});

		var error = targetEEX - currentEEX;

		var deltaQ = stepSize * jacobian.Transpose() * error;

		Arm.SetTheta1(Arm.Theta1 + (float)deltaQ[0, 0]);
		Arm.SetTheta2(Arm.Theta2 + (float)deltaQ[1, 0]);
		Arm.SetD3(Arm.D3 + (float)deltaQ[2, 0]);
		Arm.SetTheta4(Arm.Theta4 + (float)deltaQ[3, 0]);
	}



	public override void _PhysicsProcess(double delta)
	{
		Target.Position = new Vector3((float)TargetY, (float)TargetZ, (float)TargetX);
		Target.Rotation = new Vector3(0, (float)TargetTheta, 0);

		switch (IKMethod)
		{
			case IKMethod.PseudoInverse:
				PseudoInverseIK(StepSize);
				break;
			case IKMethod.DampedLeastSquaresWithGain:
				DampedLeastSquaresIKWithGain(lambda: 1.0, stepSize: StepSize);
				break;
			case IKMethod.JacobianTranspose:
				JacobianTransposeIK(stepSize: StepSize);
				break;
		}

		PlotEELocation();

	}

	private void PlotEELocation()
	{
		var eeLocation = Arm.GetCurrentEEPosition();
		MultiMesh.VisibleInstanceCount++;
		var transform = new Transform3D(Basis.Identity, eeLocation);
		MultiMesh.SetInstanceTransform(MultiMesh.VisibleInstanceCount - 1, transform);
	}
}
