package at.monol1th;

import at.monol1th.cpic.core.lattice.Lattice;

public class Main {

    public static void main(String[] args) {

	    int s = 60;
	    int n = 3;

	    Lattice testLat = new Lattice(new int[]{s,s,s},n, 0.1);

	    System.out.println(testLat.getPlaquette(1, 2, 3).getTrace().toString());

    }
}
