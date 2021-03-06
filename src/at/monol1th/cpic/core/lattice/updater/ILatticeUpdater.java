package at.monol1th.cpic.core.lattice.updater;

import at.monol1th.cpic.core.lattice.Lattice;

/**
 * Created by dmueller on 3/30/15.
 */
public interface ILatticeUpdater
{
    public void updateLattice(Lattice L, double timeStep);
}
