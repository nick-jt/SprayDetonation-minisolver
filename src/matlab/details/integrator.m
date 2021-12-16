
function [t,y,M]=integrator(U0,vars)
 
 
[Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D, lam ...
    , alpha, Length, fuel, phi, mech, q, gas...
    , satpressure, latheat, dropCv] = vars{1:end};

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
Lspan = [0 Length];
y0 = [ Tg1; Rhog1; Ug1; Td1; Ud1; Rd1; Yg1 ];

vars = {Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
    lam, alpha, Length, fuel, phi, mech, q, gas, satpressure, latheat, dropCv};

% Solving IVP
options = odeset('OutputFcn',@(t,y,flag) myoutputfunc(t,y,flag,gas),...
    'NonNegative',1:6,'RelTol',1e-5,'AbsTol',1e-8);
[t,y] = ode15s(@(t,y) statevectorfunction(t,y,vars), Lspan, y0, options);

M = zeros(length(t),1);
for i = 1:length(t)
    set(gas,'T',y(i,1),'Rho',y(i,2),'Y',y(i,7:end));
    M(i) = y(i,3)/soundspeed(gas);
end

end
