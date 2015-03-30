package at.monol1th.cpic.core.lattice;

import org.apache.commons.math3.complex.Complex;

/**
 * Created by David on 28.03.2015.
 */
public class Lattice
{
	/*
		N is the dimensionality of the lattice.
		size is an array containing the lattice sizes.
	 */
	public int dim;
	public int[] size;

	/*
		N defines the gauge group SU(N).
	 */
	public int N;

	/*
		Spatial distance between two lattice neighbouring points.
	 */
	public double gridSpacing;

    /*
        Yang-Mills coupling constant
     */
    public double couplingConstant;

	/*
		Array of cells. A cell sits at each lattice point and contains the gauge links and momenta.
	 */
	public Cell[] cells;

	/*
		Total number of cells in the lattice.
	 */
	private int numberOfCells;

	/*
		Array which translates the cell index to its spatial position in the lattice.
	 */
	private int[][] indexToPositionArray;

	/*
		Primary constructor of the Lattice class.
	 */
	public Lattice(int[] size, int N, double gridSpacing, double couplingConstant)
	{
		this.size = size;
		this.gridSpacing = gridSpacing;
        this.couplingConstant = couplingConstant;
		this.dim = size.length;

		this.N = N;

		/*
			Compute total number of cells. Duh.
		 */
		this.numberOfCells = 1;
		for(int i = 0; i < size.length; i++)
		{
			numberOfCells *= size[i];
		}

		/*
			Initialize lattice of cells.
		 */
		generateIndexToPositionArray();
		cells = new Cell[numberOfCells];
		for(int i = 0; i < numberOfCells; i++)
		{
			int[] pos = indexToPositionArray[i];


			/*
				Setup array of neighbouring cell indices.
			 */
			int[] neighbors = new int[2*dim];
			for(int j = 0; j < dim; j++)
			{
				int[] shifted_pos = pos.clone();
				shifted_pos[j] += 1;
				neighbors[2*j] = getCellIndex(shifted_pos);
				shifted_pos[j] -= 2;
				neighbors[2*j + 1] = getCellIndex(shifted_pos);
			}

			cells[i] = new Cell(N, dim, indexToPositionArray[i], neighbors);
		}

	}

    /*
        Electric field methods
     */

    public ComplexMatrix getElectricField(Cell cell, int direction)
    {
        ComplexMatrix c1 = ComplexMatrix.multiply(getGaugeLink(cell, direction).getHermitianConjugate(), getGaugeLinkDerivative(cell, direction));
        ComplexMatrix c2 = ComplexMatrix.multiply(getGaugeLinkDerivative(cell, direction), getGaugeLink(cell, direction).getHermitianConjugate());

        Complex factor = new Complex(0.0, 1.0/(2.0 * couplingConstant * gridSpacing));

        return ComplexMatrix.multiply(ComplexMatrix.add(c1, c2), factor);
    }

    public ComplexMatrix getElectricField(int pos[], int direction)
    {
        return getElectricField(getCell(pos), direction);
    }

    public ComplexMatrix getElectricField(int index, int direction)
    {
        return getElectricField(getCell(index), direction);
    }


	/*
		Plaquette methods
	 */

	public ComplexMatrix getPlaquette(Cell c, int i, int j)
	{
		Cell c2 = getNeighbouringCell(c, i);
		Cell c3 = getNeighbouringCell(c2, j);
		Cell c4 = getNeighbouringCell(c3, -i);

		ComplexMatrix l1 = getGaugeLink(c, i);
		ComplexMatrix l2 = getGaugeLink(c2, j);
		ComplexMatrix l3 = getGaugeLink(c3, -i);
		ComplexMatrix l4 = getGaugeLink(c4, -j);

		return ComplexMatrix.multiply(new ComplexMatrix[] {l4,l3,l2,l1});
	}

	public ComplexMatrix getPlaquette(int pos[], int i, int j)
	{
		return getPlaquette(getCell(pos), i, j);
	}

	public ComplexMatrix getPlaquette(int index, int i, int j)
	{
		return getPlaquette(getCell(index), i, j);
	}

    /*
        getStaple methods
     */

    public ComplexMatrix getStaple(Cell c, int i, int j)
    {
        Cell c2 = getNeighbouringCell(c, i);
        Cell c3 = getNeighbouringCell(c2, j);
        Cell c4 = getNeighbouringCell(c3, -i);

        ComplexMatrix l2 = getGaugeLink(c2, j);
        ComplexMatrix l3 = getGaugeLink(c3, -i);
        ComplexMatrix l4 = getGaugeLink(c4, -j);

        return ComplexMatrix.multiply(new ComplexMatrix[] {l4,l3,l2});
    }

    public ComplexMatrix getStaple(int pos[], int i, int j)
    {
        return getStaple(getCell(pos), i, j);
    }

    public ComplexMatrix getStaple(int index, int i, int j)
    {
        return getStaple(getCell(index), i, j);
    }

	/*
		getGaugeLink methods
	 */

	public ComplexMatrix getGaugeLink(Cell c, int direction)
	{
		if(direction > 0)
		{
			return c.getGaugeLink(direction);
		}
		else
		{
			return getNeighbouringCell(c, direction).getGaugeLink(-direction).getHermitianConjugate();
		}
	}

	public ComplexMatrix getGaugeLink(int[] pos, int direction)
	{
		return getGaugeLink(getCell(pos), direction);
	}

	public ComplexMatrix getGaugeLink(int index, int direction)
	{
		return getGaugeLink(getCell(index), direction);
	}

    /*
        getGaugeLinkDerivative methods.
            Returns \dot{U} from momentum \Pi
     */

    public ComplexMatrix getGaugeLinkDerivative(Cell cell, int direction)
    {
        if(direction > 0)
        {
            return ComplexMatrix.multiply(cell.getMomentum(direction), Math.pow(couplingConstant, 2.0) / Math.pow(gridSpacing, dim-2)).getConjugate();
        }
        else
        {
            return ComplexMatrix.multiply(getNeighbouringCell(cell, direction).getMomentum(-direction), Math.pow(couplingConstant, 2.0) / Math.pow(gridSpacing, dim-2)).getConjugate();
        }
    }

    public ComplexMatrix getGaugeLinkDerivative(int[] pos, int direction)
    {
        return getGaugeLinkDerivative(getCell(pos), direction);
    }

    public ComplexMatrix getGaugeLinkDerivative(int index, int direction)
    {
        return getGaugeLinkDerivative(getCell(index), direction);
    }

	/*
		getCell methods
	 */

	public Cell getCell(int[] pos)
	{
		return cells[getCellIndex(pos)];
	}

	public Cell getCell(int index)
	{
		return cells[index];
	}

	/*
		getNeighbouringCell methods
	 */

	public Cell getNeighbouringCell(int[] pos, int direction)
	{
		return getNeighbouringCell(getCell(pos), direction);
	}

	public Cell getNeighbouringCell(int index, int direction)
	{
		return getNeighbouringCell(getCell(index), direction);
	}

	public Cell getNeighbouringCell(Cell c, int direction)
	{
		return getCell(c.getNeighbourIndex(direction));
	}

	/*
		Methods related to constructing index arrays and computing cell indices and positions.
	 */

	private void generateIndexToPositionArray()
	{
		this.indexToPositionArray = new int[numberOfCells][dim];
		 for(int i = 0 ; i < numberOfCells; i++)
		 {
			 indexToPositionArray[i] = getPosition(i);
		 }

	}

	private int getCellIndex(int[] pos)
	{
		for(int i = 0; i < dim; i++)
		{
			pos[i] = (pos[i] % size[i] + size[i]) % size[i];
		}

		int cellPos = pos[0];

		for(int i = 1; i < pos.length; i++)
		{
			cellPos *= size[i];
			cellPos += pos[i];
		}
		return  cellPos;
	}

	private int[] getPosition(int cellIndex)
	{
		int[] pos = new int[dim];

		for(int i = dim-1; i >= 0; i--)
		{
			pos[i] = cellIndex % size[i];
			cellIndex -= pos[i];
			cellIndex /= size[i];
		}

		return pos;
	}

}
