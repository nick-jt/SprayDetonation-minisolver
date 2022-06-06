import sys
import cantera as ct
import pyjacob
import numpy as np
import numpy.linalg as ln


if __name__ == '__main__':

        gas = ct.Solution('Dodecane.cti')
        #for i,arg in enumerate(sys.argv):
        #       print(i,arg)
        n2_ind = gas.species_index('N2')
        specs = gas.species()[:]
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                species=specs[:n2_ind] + specs[n2_ind + 1:] + [specs[n2_ind]],
                reactions=gas.reactions())
        gas.TPY = float(sys.argv[1]),float(sys.argv[2]),[float(v) for v in sys.argv[3:]]

        #setup the state vector
        y = np.zeros(gas.n_species)
        y[0] = gas.T
        y[1:] = gas.Y[:-1]

        #create a dydt vector
        #dydt = np.zeros_like(y)
        #       pyjacob.py_dydt(0, P, y, dydt)

        #create a jacobian vector
        jac = np.zeros(gas.n_species * gas.n_species)

        #evaluate the Jacobian
        pyjacob.py_eval_jacobian(0, gas.P, y, jac)

        jac = np.reshape(jac,(int(np.sqrt(len(jac))), int(np.sqrt(len(jac)))))
        eigs = ln.eigvals(jac)
        #for i,eig in enumerate(eigs):
        #       print(gas.species_name(i),1/eig)

        #print("inverse eig timescale ",min(abs(1/eigs)))
        scales = abs(1/eigs)
        COscale = scales[gas.species_index('CO')+1]
        OHscale = scales[gas.species_index('OH')+1]
        Fuelscale = scales[gas.species_index('NC12H26')+1]
        Tscale = scales[0]
        eigenscale = min(abs(1/eigs)) # inverse eigenvalue timescale
        ReacTscales = [eigenscale, COscale, OHscale, Fuelscale, Tscale]
