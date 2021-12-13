function [dydx,extras] = statevectorfunction(t,y,vars)
%% Basic variables


if (y(1)<0) y(1)=1;
elseif (isnan(y(1))) y(1)=1;
end

if (y(2)<0) y(2)=1;
elseif (isnan(y(2))) y(2)=1;
end

% Unpacking evolution values
Tg = y(1); rhog = y(2); ug = y(3); Td = y(4); ud = y(5); rd = y(6);
Yg = y(7:end);
 
[~, ~, Cdw, Chw, Rd0, l_char, Pr, Le, Tw, Cvd, rhod, nu0, D,...
    lam, ~, ~, fuel, ~, ~, ~, gas] = vars{1:end};
 
% Get vapor enthalpy at droplet temperature
fuel_index = speciesIndex(gas,fuel);
w_k = molecularWeights(gas);
if (rd>1e-2*Rd0)
    set(gas,'T',Td,'Rho',rhog,'Y',Yg);
    enth = enthalpies_RT(gas)*gasconstant*Td/w_k(fuel_index);
    hgf_Td = enth(fuel_index);
else
    hgf_Td = 0;
end

% Reset gas state
set(gas,'T',Tg,'Rho',rhog,'Y',Yg);
 
% Gas Parameters
Cpg     = cp_mass(gas);
w       = meanMolecularWeight(gas);
hg_k    = enthalpies_RT(gas)*gasconstant*Tg./w_k;
hgf_Tg  = hg_k(fuel_index);
omega   = netProdRates(gas);
gamma   = Cpg/cv_mass(gas);
grs     = gamma/(gamma-1);
wf      = w_k(fuel_index);
nd      = nu0/ud;
Re      = rhog * abs(ud-ug) * 2*rd/viscosity(gas);
conv    = ( 1 + 0.276*Re^0.5*Pr^0.5 );
wnf     = sum(w_k(2:end).*Yg(2:end))/sum(Yg(2:end)); % TODO
M       = ug/soundspeed(gas);

%% Droplet Empirical Equations
if (rd>1e-2*Rd0)
    
    % Latent Heat of Vaporization
    A = 53.66; B = 0.2831; Tc = 540.2;
    Tr = Td/Tc; 
    L = A*exp(-B*Tr) * (1-Tr)^B;
    L = L / wf * 1000 * 1000; % J/kg(fuel) TODO

    % Vaporization Rate (Antoine Equation)
    % pfs2 = satPressure(gas,Tg); 
    A = 4.02832; B = 1268.636; C = -56.199;
    pfs = 10^(A - B/(Td + C)) * 100000;
    Xfs = pfs/pressure(gas);
    Yfs = Xfs*wf/(Xfs*wf+(1-Xfs)*wnf);
    By = (Yfs-Yg(fuel_index))/(1-Yfs);
    Bh = Cpg*(Tg-Td)/L;
    CDd = 22*Re^-1 * conv;

    fd = nd*CDd*4*pi*rd^2*rhog*abs(ud-ug)*(ud-ug)/2;
    mdotv = nd*4*pi*rd*lam/(Le*Cpg) * log(1+By)*conv;
    qd = nd*4*pi*rd*lam/Cpg*log(1+Bh)*conv*L;
    
    % Droplet Radius
    drddx = -mdotv/(rhod*4*pi*rd^2*nu0);

    % Droplet Velocity
    duddx = -fd/(rhod*nu0*4/3*pi*rd^3);

    % Droplet Temperature
    dTddx = (qd-mdotv*L)/(rhod*nu0*4/3*pi*rd^3*Cvd);
    
else
    [fd,mdotv,qd,drddx,duddx,dTddx] = deal(0);
end 
    
%% Wall Losses
fw = Cdw/l_char * rhog * abs(D-ug)*(D-ug)/2;
qw = Chw/l_char * rhog*abs(D-ug)*Cpg*(Tg-Tw);

%% Evolution parameters

% Droplet Vaporization
enthalpy_from_vaporization = hgf_Tg-hgf_Td;

% Gas Velocity
ExternalLosses  = fw*(grs*ug-D) + qw + fd*(grs*ug-ud)+qd;
Thermicity      = sum((w./w_k - hg_k./(Cpg*Tg)).*omega.*w_k);
Contr_From_Dropl= mdotv * (grs*ug*(ud-ug)-Cpg*Tg*w/wf...
    -(ud^2-ug^2)/2 + enthalpy_from_vaporization);
S = ExternalLosses-Cpg*Tg*Thermicity+Contr_From_Dropl;

dugdx = (gamma-1)*M^2*S/((M^2-1)*rhog*ug^2);

% Gas Temperature
dTgdx = -ug/Cpg*dugdx + 1/(rhog*ug*Cpg) * ...
    (fw*D + fd*ud - sum(hg_k.*omega.*w_k) - qw - qd ...
    + mdotv*((ud^2-ug^2)/2 - enthalpy_from_vaporization));

% Gas Density
drhogdx = -rhog/ug*dugdx + mdotv/ug;

% Species evolution
Ydk = zeros(nSpecies(gas),1);
Ydk(fuel_index) = 1;
dYgdx = 1/(rhog*ug) * (omega.*w_k+mdotv*(Ydk-Yg));

dydx = [dTgdx; drhogdx; dugdx; dTddx; duddx; drddx; dYgdx];
extras = [M; Thermicity];
end
