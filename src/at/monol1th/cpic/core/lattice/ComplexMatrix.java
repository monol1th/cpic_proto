package at.monol1th.cpic.core.lattice;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

/**
 * Created by David on 28.03.2015.
 */
public class ComplexMatrix {

	public RealMatrix real;
	public RealMatrix imag;

	public int N;

	public ComplexMatrix(int N, boolean unit)
	{
		this.N = N;
		this.real = new Array2DRowRealMatrix(N, N);
		this.imag = new Array2DRowRealMatrix(N, N);
		if(unit)
			real = MatrixUtils.createRealIdentityMatrix(N);
	}

	public ComplexMatrix(RealMatrix real, RealMatrix imag)
	{
		this.N = real.getRowDimension();
		this.real = real.copy();
		this.imag = imag.copy();
	}

	public ComplexMatrix(ComplexMatrix comp)
	{
		set(comp);
	}

	public void set(ComplexMatrix comp)
	{
		this.N = comp.N;
		this.real = comp.real.copy();
		this.imag = comp.imag.copy();
	}

	public static ComplexMatrix add(ComplexMatrix a, ComplexMatrix b)
	{
		ComplexMatrix c = a.copy();
		c.real.add(b.real);
		c.imag.add(b.imag);
		return c;
	}

	public static ComplexMatrix add(ComplexMatrix[] matrices)
	{
		ComplexMatrix c = matrices[0].copy();
		for(int i = 1; i < matrices.length; i++)
		{
			c = ComplexMatrix.add(c, matrices[i]);
		}
		return c;
	}

	public static ComplexMatrix multiply(ComplexMatrix a, ComplexMatrix b)
	{
		ComplexMatrix c1 = a.copy();
		ComplexMatrix c2 = a.copy();

		c1.real.multiply(b.real);
		c1.real.subtract(c2.imag.multiply(b.imag));

		c1.imag.multiply(b.real);
		c1.imag.add(c2.real.multiply(b.imag));

		return new ComplexMatrix(c1);
	}

	public static ComplexMatrix multiply(ComplexMatrix[] matrices)
	{
		ComplexMatrix c = matrices[0].copy();
		for(int i = 1; i < matrices.length; i++)
		{
			c = ComplexMatrix.multiply(c, matrices[i]);
		}
		return c;
	}

	public Complex getTrace()
	{
		return new Complex(real.getTrace(), imag.getTrace());
	}

	public ComplexMatrix getHermitianConjugate()
	{
		// TODO: asdaosid
		return new ComplexMatrix(real.copy().transpose(), imag.copy().scalarMultiply(-1.0).transpose());
	}

	public ComplexMatrix copy()
	{
		return new ComplexMatrix(this);
	}
}
