# SprayDetonation-minisolver
A multi-phase ZND detonation solver. 

## Testing
To test, start MATLAB and run the following code.
```
cd SprayDetonation-minisolver/src/matlab
[x,y,M] = SprayDet_Matlab;
```

`x` is a vector of distances from the shock. `M` is a mach number vector. `y` is a 2-d matrix of size `[6 + num_species]x[num_timesteps]` containing the following parameters.

1. gas temperature
2. gas density
3. gas velocity
4. droplet temperature
5. droplet velocity
6. droplet radius
7. and beyond. Species mass fractions.

After running the code and loading the solution into the MATLAB variable space. The following code will plot gas temperature with distance from shock.

```
plot(x,y(:,1))
```

## Case conditions
CaseSetup.m contains most of the relevant parameters for your detonation. Change equivalence ratio, droplet loading, droplet diameter, external loss coefficients, initial conditions, etc...

The `funcflag` in CaseSetup.m defines the solution method for the detonation. SprayDet_Dodecane.m then will begin the bracketing method iteration by calling the respective matlab file. bracketMethodCDW.m file iterates to determine the drag coefficient for a given detonation velocity (may be used to capture both stable and unstable solutions of the detonation velocity vs loss curve). You may also use getSSvelocity.m to iterate for the detonation steady state velocity given case conditions, including heat and friction wall loss coefficients.
