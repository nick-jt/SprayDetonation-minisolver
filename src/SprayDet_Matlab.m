function SprayDet_Matlab

addpath(genpath('/home/njt14003/Cantera/cantera-build/lib/cantera/matlab/toolbox'))
addpath('details');
addpath('mechanisms');

% Most common parameters
phi     = 0.7;
Cdw     = 0.00;
Chw     = 0.00;
Rd0     = 5e-6;
alpha   = 1;
T0      = 298.15;
P0      = 1e5;
fuel    = 'C7H16';     % Change parameters if this changes
mech    = 'Heptane.cti';


% Parameters
lchar   = 3.81*0.01/4;
Pr      = 1;
Le      = 1;
Tw      = 298.15;
Cvd     = 2236;
rhod    = 680;
Length  = 0.2;
nu0     = 0;
lam     = 0;

gas     = Solution(mech);
C_count = nAtoms(gas,fuel,'C');
H_count = nAtoms(gas,fuel,'H');
a       = C_count + 0.25 * H_count;			
q       = fuel + ":" + string(phi*(1-alpha)) ...
	+ ", O2:" + string(a) ...
	+ ", N2:" + string(a*3.76);
    
vars = {T0 P0 Cdw Chw Rd0 lchar Pr Le...
Tw Cvd rhod nu0 U0 lam alpha Length...
fuel phi mech char(q) gas};

[t,y,M] = integrator(U0,vars);

rmpath('details');
rmpath('mechanisms');

end
