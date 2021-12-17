% Use this file to initialize the parameters for your case

function [Vtests,SS_Cdws] = SprayDet_Matlab_JP10
clear

addpath('details');

% Most common parameters
phi     = 1.0;
Cdw     = 0.00;
Chw     = 0.00;
Rd0     = 5e-6;
alpha   = 1;      % initial droplet loading [kg liquid fuel/ kg total fuel]
T0      = 298.15; % initial T [K]
P0      = 1e5;    % initial P [Pa]
fuel    = 'C10H16';     % Change parameters if this changes
mech    = 'JP10.cti';

% Droplet related parameters
satpressure = @(Td) exp(8.173667-2783.0/(Td-117.037))*101325; % Antoine Equation
latheat = @(Td,wf) 541911.764706;
dropCv = @(Td,w) gasconstant/w*(3.3218+0.07975*Td^0.85+27.6975*(1470/Td)^2 ...
    *exp(1470/Td)/(exp(1470/Td)-1)^2);

% Parameters
lchar   = 3.81*0.01/4;
Pr      = 1;
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
q       = fuel + ":" + string(phi*(1-alpha)) ...
        + ", O2:" + string(a) ...
        + ", N2:" + string(a*3.76);

U0 = 1200;
vars = {T0 P0 Cdw Chw Rd0 lchar Pr Le...
	Tw Cvd rhod nu0 U0 lam alpha Length...
	fuel phi mech char(q) gas satpressure latheat dropCv};

printcase(vars);

% [t,y,M] = integrator(U0,vars);
% Vtests = linspace(1200,1800,20);
% SS_Cdws = zeros(1,length(Vtests));
% for i=1:length(Vtests)
%     fprintf("V=%f\n",Vtests(i));
    [SS_Cdws(i),~,~,~] = bracketMethodCDW(U0,0,0.1,vars,1);
% end
% [SS_Velocity,t,y,M] = getSSvelocity(1700,1900,vars,1);


% [Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
%     lam, alpha, Length, fuel, phi, mech, q, gas, satpressure, latheat] = vars{1:end};
% set(gas, 'T', Tg0, 'P', Pg0, 'X', char(q));
% rho0 = density(gas);
% if (Rd0>0 && alpha>0)
%     nd = getnumden(gas,alpha,phi,rhod,Rd0,fuel,rho0);
% else
%     nd = 0;
% end
% set(gas, 'T', Tg0, 'P', Pg0, 'X', char(q));
% gas = postshockstate( U0, Pg0, Tg0, char(q), mech );
% lam = thermalConductivity(gas);
% Tg1 = temperature(gas);
% Rhog1 = density(gas);
% Ug1 = rho0*U0/Rhog1;
% Td1 = Tg0;
% Ud1 = U0;
% Rd1 = Rd0;
% Yg1 = massFractions(gas);
% nu0 = nd*Ud1;
% vars = {Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
%     lam, alpha, Length, fuel, phi, mech, q, gas, satpressure, latheat, dropCv};
% extras = zeros(length(t),6);
% for i = 1:length(t)
%     extras(i,:) = getQsrcterms(t,y(i,:),vars);
% end
% 
% if (size(t)==1)
%     error("size(t)=1")
% end

rmpath('details');

end
