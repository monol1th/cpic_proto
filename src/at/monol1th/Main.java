package at.monol1th;

import at.monol1th.cpic.core.lattice.Lattice;

public class Main {

    public static void main(String[] args) {

	    int s = 10;
	    int n = 2;

	    Lattice testLat = new Lattice(new int[]{s,s,s},n, 0.1, 0.1);

	    System.out.println(testLat.getElectricField(1, 2).getTrace().toString());
        System.out.println(testLat.getGaugeLinkDerivative(10, 2).getTrace().toString());

    }
}
