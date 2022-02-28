function [dydx,extras] = statevectorfunction(t,y,vars)
%% Basic variables


% Unpacking evolution values
Tg = y(1); rhog = y(2); ug = y(3); Td = y(4); ud = y(5); rd = y(6);
Yg = y(7:end);
 
[~, ~, Cdw, Chw, Rd0, l_char, Pr, Le, Tw, ~, rhod, nu0, D,...
    lam, ~, ~, fuel, ~, ~, ~, gas, satpressure, latheat, dropCv] = vars{1:end};

% Use gas state to get vapor enthalpy at droplet temperature and mole fractions without fuel
fuel_index = speciesIndex(gas,fuel);
w_k = molecularWeights(gas);
if (rd>1e-2*Rd0 && nu0>0)
    set(gas,'T',Td,'Rho',rhog,'Y',Yg);
    enth = enthalpies_RT(gas)*gasconstant*Td/w_k(fuel_index);
    hgf_Td = enth(fuel_index);
else
    hgf_Td = 0;
end

% Mean molecular weight with no fuel
fidx = speciesIndex(gas,fuel);
if (massFraction(gas,fidx)>0)
	Yg_nofuel = Yg;
	Yg_nofuel(fidx) = 0;
	set(gas,'T',Tg,'Rho',rhog,'Y',Yg);
	wnf = meanMolecularWeight(gas);
else
	wnf = meanMolecularWeight(gas);
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
mu	= viscosity(gas);
Re      = rhog * abs(ud-ug) * 2*rd/mu;
conv    = ( 1 + 0.276*Re^0.5*Pr^0.5 );
c	= soundspeed(gas);
M       = ug/c;

%% Droplet Empirical Equations
if (rd>1e-2*Rd0 && nu0>0)
    Tb = 489;

    Lv = latheat(Td,wf);
    Cvd = dropCv(Td,w);
    
    % film temp properties
    Cpf = cp_mass(gas);
    kf = thermalConductivity(gas);
    
    % Mass transfer
    R = gasconstant();
    P = pressure(gas);
    if ~exist('taum','var')
        Xeq = 101325/P * exp(Lv/(R/wf)*(1/Tb-1/Td));
        Xeq = min(1-1e-10,Xeq);
        Yeq = Xeq / (Xeq + (1-Xeq)*w/wf);
        By = (Yeq - massFraction(gas,fuel)) / (1 - Yeq);
        Sh = 2+0.552*Re^0.5*Pr^(1/3);
        global taum;
        taum =  (4*rd^2)*rhod/(6*Sh*kf/(Le*Cpf)*log(1+By));
    end
    beta = rhod*Cpf*(rd*2)^2/(12*kf*taum); % TODO get film temp
    Lk = kf/(Le*Cpf)*sqrt(2*pi*Td*R/wf)/P;
    Xeq = 101325/P * exp(Lv/(R/wf)*(1/Tb-1/Td));
    Xeq = min(1-1e-10,Xeq);
    Xneq = Xeq - Lk/rd*beta;
    fprintf("%f %f %f %f\n",Xeq,Lk,beta,rd);
    Yneq = Xneq / (Xneq + (1-Xneq)*w/wf);
    By = (Yneq - massFraction(gas,fuel)) / (1 - Yneq);
    taum1 =  (4*rd^2)*rhod/(6*Sh*kf/(Le*Cpf)*log(1+By));
    mdotv = nd*4*pi*rd*lam/(Le*Cpg) * log(1+By)*conv;
    taum = taum1;
     
    CDd = dragcoefficient(y,mu,c);
    fd = nd*CDd*pi*rd^2*rhog*abs(ud-ug)*(ud-ug)/2;
    
    Bh = Cpg*(Tg-Td)/Lv;
    qd = nd*4*pi*rd*lam/Cpg*log(1+Bh)*conv*Lv;
    
    % Droplet Radius
    drddx = -mdotv/(rhod*4*pi*rd^2*nu0);
 
    % Droplet Velocity
    duddx = -fd/(rhod*nu0*4/3*pi*rd^3);

    % Droplet Temperature
    dTddx = (qd-mdotv*Lv)/(rhod*nu0*4/3*pi*rd^3*Cvd);
    
else
    [fd,mdotv,qd,drddx,duddx,dTddx] = deal(0);
end 
    
%% Wall Losses
fw = Cdw/l_char * rhog * abs(D-ug) * (D-ug)/2;
qw = Chw/l_char * rhog * abs(D-ug) * Cpg * (Tg-Tw);

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
