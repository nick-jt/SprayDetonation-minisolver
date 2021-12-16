% Use this file to initialize the parameters for your case

function [t,y,M,extras] = SprayDet_Matlab_JP10
clear

addpath('details');

% Most common parameters
phi     = 1.0;
Cdw     = 0.00;
Chw     = 0.00;
Rd0     = 0;
alpha   = 0;
T0      = 298.15;
P0      = 1e5;
fuel    = 'C10H16';     % Change parameters if this changes
mech    = 'JP10.cti';

% Droplet related parameters
satpressure = @(Td) exp(8.173667-2783.0/(Td-117.037))*101325; % Antoine Equation
latheat = @(Td,wf) 541911.764706;
dropCv = @(Td,w) gasconstant/w*(3.3218+0.07975*Td^0.85+27.6975*(1470/Td)^2 ...
    *exp(1470/Td)/(exp(1470/Td)-1)^2);

% Parameters
lchar   = 3.81*0.01/4;
Pr      = 0.72;
Le      = 1;
Tw      = T0;
Cvd     = 0; % replaced later
rhod    = 932;
Length  = 0.2;
nu0     = 0; % replaced later
lam     = 0; % replaced later

gas     = Solution(mech);
C_count = nAtoms(gas,fuel,'C');
H_count = nAtoms(gas,fuel,'H');
a       = C_count + 0.25 * H_count;			
% q       = fuel + ":" + string(phi*(1-alpha)) ...
% 	+ ", O2:" + string(a) ...
% 	+ ", N2:" + string(a*3.76);
    
q       = fuel + ":" + string(phi*(1-alpha)) ...
        + ", O2:" + string(a);

U0 = 1770;
vars = {T0 P0 Cdw Chw Rd0 lchar Pr Le...
	Tw Cvd rhod nu0 U0 lam alpha Length...
	fuel phi mech char(q) gas satpressure latheat dropCv};

printcase(vars);

% [t,y,M] = integrator(U0,vars);
% [SS_Cdw,t,y,M] = bracketMethodCDW(V,CDWlow,CDWhigh,vars,fileout)
[SS_Velocity,t,y,M] = getSSvelocity(1900,2500,vars,1);


[Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
    lam, alpha, Length, fuel, phi, mech, q, gas, satpressure, latheat] = vars{1:end};

% Initialize the gas
set(gas, 'T', Tg0, 'P', Pg0, 'X', char(q));
rho0 = density(gas);

% Calculate number density (modifies gas state)
if (Rd0>0 && alpha>0)
    nd = getnumden(gas,alpha,phi,rhod,Rd0,fuel,rho0);
else
    nd = 0;
end

% Initialize the gas
set(gas, 'T', Tg0, 'P', Pg0, 'X', char(q));
 
% Calculate Postshock Conditions
gas = postshockstate( U0, Pg0, Tg0, char(q), mech );
lam = thermalConductivity(gas);

% IVP initial conditions
Tg1 = temperature(gas);
Rhog1 = density(gas);
Ug1 = rho0*U0/Rhog1;
Td1 = Tg0;
Ud1 = U0;
Rd1 = Rd0;
Yg1 = massFractions(gas);

% Initializing
nu0 = nd*Ud1;
vars = {Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
    lam, alpha, Length, fuel, phi, mech, q, gas, satpressure, latheat, dropCv};


extras = zeros(length(t),6);
for i = 1:length(t)
    extras(i,:) = getQsrcterms(t,y(i,:),vars);
end

if (size(t)==1)
    error("size(t)=1")
end

rmpath('details');

end
