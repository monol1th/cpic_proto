package at.monol1th.cpic.core.lattice;
/**
 * Created by David on 28.03.2015.
 */
public class Cell {
	public int N;
	public int d;
	public int[] x;

	public ComplexMatrix[] links;
	public ComplexMatrix[] momenta;

	public int[] neighbouringCells;

	public Cell(int N, int d, int[] x, int[] neighbours)
	{
		this.N = N;
		this.d = d;
		this.x = x;

		links = new ComplexMatrix[d];
		for(int i = 0; i < d; i++)
			links[i] = new ComplexMatrix(N, true);

		momenta = new ComplexMatrix[d];
		for(int i = 0; i < d; i++)
			momenta[i] = new ComplexMatrix(N, true);

		this.neighbouringCells = neighbours;
	}

	public ComplexMatrix getGaugeLink(int d)
	{
		return links[d-1];
	}

	public ComplexMatrix getMomentum(int d)
	{
		return momenta[d-1];
	}

	public int getNeighbourIndex(int direction)
	{
		int i = Math.abs(direction) - 1;
		int j = (1 - Integer.signum(direction)) / 2;

		return neighbouringCells[2*i+j];
	}
}
