function [dydx,extras] = statevectorfunction(t,y,vars)

% Unpacking evolution values
Tg = y(1); rhog = y(2); ug = y(3); Td = y(4); ud = y(5); rd = y(6);
Yg = y(7:end);
 
[~, ~, Cdw, Chw, Rd0, l_char, Pr, Le, Tw, ~, rhod, nu0, D,...
    lam, ~, ~, fuel, ~, ~, ~, gas, satpressure, latheat, dropCv] = vars{1:end};


% set gas state
try
	set(gas,'T',Tg,'Rho',rhog,'Y',Yg);
catch ERROR
	fprintf(['ERROR: ' ERROR.message]);
	fprintf("	ending run and outputting data");
	dydx=zeros(length(y),1);
	extras=zeros(5,1);;
	return
end

% Gas Parameters
fuel_index = speciesIndex(gas,fuel);
w_k 	= molecularWeights(gas);
P 	= pressure(gas);
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
c	= soundspeed(gas);
M       = ug/c;

%% Using gas state, will never be needed again

% Use gas state to get vapor enthalpy at droplet temperature and mole fractions without fuel
if (rd>1e-2*Rd0 && nu0>0) % TODO use film temperature?
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
	set(gas,'P',P,'Rho',rhog,'Y',Yg_nofuel);
	wnf = meanMolecularWeight(gas);
else
	wnf = meanMolecularWeight(gas);
end

global liquid_vaporized_flag;

%% Droplet Empirical Equations
if (rd>1e-3*Rd0 && nu0>0)

    liquid_vaporized_flag = 0;

    Tb = 489;

    Lv = latheat(Td,wf);
    Cvd = dropCv(Td,w);
    
    % film temp properties TODO
    Twb = 137*(Tb/373.15)^0.68 * log10(Tg) - 45;
    set(gas,'T',Twb,'Rho',rhog,'Y',Yg);
    Cpf = cp_mass(gas);
    kf = thermalConductivity(gas);
    muf = viscosity(gas);
    rhof = density(gas);
    Ref   = rhof * abs(ud-ug) * 2*rd/muf;
    Prf   = muf*Cpf/kf;
    conv = ( 1 + 0.276*Ref^0.5*Prf^(1/3) ); % Sc=Pr for Le=1
    
    % Mass transfer
    R = gasconstant();
    if ~exist('taum','var') % For our first iteration
        Xeq = 101325/P * exp(Lv/(R/wf)*(1/Tb-1/Td));
        Xeq = min(1-1e-3,Xeq); %TODO need this?
        Yeq = Xeq / (Xeq + (1-Xeq)*w/wf);
        By = (Yeq - massFraction(gas,fuel)) / (1 - Yeq);
        Sh = 2*conv;
        global taum;
        taum =  (4*rd^2)*rhod/(6*Sh*kf/(Le*Cpf)*log(1+By));
    end
    beta = rhod*Cpf*(rd*2)^2/(12*kf*taum); 
    Lk = kf/(Le*Cpf)*sqrt(2*pi*Td*R/wf)/P;
    Xeq = 101325/P * exp(Lv/(R/wf)*(1/Tb-1/Td)); 
    Xeq = min(1-1e-3,Xeq); % TODO Need this?
    Xneq = Xeq - Lk/rd*beta;
    Yneq = Xneq / (Xneq + (1-Xneq)*w/wf);
    Yneq = min(1-1e-3,Yneq); % TODO do we need this?
    By = (Yneq - massFraction(gas,fuel)) / (1 - Yneq);
    taum1 =  (4*rd^2)*rhod/(6*Sh*kf/(Le*Cpf)*log(1+By));
    mdotv = nd*4*pi*rd*kf/(Le*Cpf) * log(1+By)*conv;
    taum = taum1;
    
    % Heat Transfer
    Nu = conv*2;
    taut = 2*taum*(exp(beta)-1)/Nu * Cvd/Cpf;
    dTddx = 1/(taut*ud)*(Tg-Td-2*Lv/Cpf*(exp(beta)-1)/Nu);
    qd = dTddx*(rhod*nd*4/3*pi*rd^3*ud*Cvd) + mdotv*Lv;
     
    % Momentum Transfer
    CDd = dragcoefficient(y,mu,c);
    fd = nd*CDd*pi*rd^2*rhog*abs(ud-ug)*(ud-ug)/2;
    
    % Droplet Radius
    drddx = -mdotv/(rhod*4*pi*rd^2*nu0);
 
    % Droplet Velocity
    duddx = -fd/(rhod*nu0*4/3*pi*rd^3);

    % Droplet Temperature
    %dTddx = (qd-mdotv*Lv)/(rhod*nu0*4/3*pi*rd^3*Cvd);
    
else
    [fd,mdotv,qd,drddx,duddx,dTddx] = deal(0);
    if ( liquid_vaporized_flag == 0)
	% Dump remaining fuel into gas mixture
	
	addedmass = zeros(1,nSpecies(gas));
	addedmass(fuel_index) = 4/3*pi*rd^3*nd*rhod/rhog;
	for i=1:nSpecies(gas)
		Yf_new = (Yg(i)+addedmass(i))/(1+addedmass(fuel_index));
	end
    end


    liquid_vaporized_flag = 1;
end 

% Resetting gas state
set(gas,'P',P,'Rho',rhog,'Y',Yg);
    
%% Wall Losses
fw = Cdw/l_char * rhog * abs(D-ug) * (D-ug)/2;
qw = Chw/l_char * rhog * abs(D-ug) * Cpg * (Tg-Tw);

%%%%%%%%%%%%%%%%%%%%%% Gas evolution

% Droplet Vaporization
enthalpy_from_vaporization = hgf_Tg-hgf_Td;

% Gas Velocity
ExternalLosses  = fw*(grs*ug-D) + qw + fd*(grs*ug-ud)+qd;
gas_thermicity      = sum((w./w_k - hg_k./(Cpg*Tg)).*omega.*w_k);
Contr_From_Dropl= mdotv * (grs*ug*(ud-ug)-Cpg*Tg*w/wf...
    -(ud^2-ug^2)/2 + enthalpy_from_vaporization);
S = ExternalLosses-Cpg*Tg*gas_thermicity+Contr_From_Dropl;
drop_thermicity = -(gamma-1)*S/(rhog*c^2);
sonicity = 1-M^2;

dugdx = drop_thermicity/sonicity;

% Gas Temperature
dTgdx = -ug/Cpg*dugdx + 1/(rhog*ug*Cpg) * ...
    (fw*D + fd*ud - sum(hg_k.*omega.*w_k) - qw - qd ...
    + mdotv*((ud^2-ug^2)/2 - enthalpy_from_vaporization));
%dPdx = -rhog*ug*dugdx+fw+fd+mdotv*(ud-ug);

% Gas Density
drhogdx = -rhog/ug*dugdx + mdotv/ug;

% Species evolution
Ydk = zeros(nSpecies(gas),1);
Ydk(fuel_index) = 1;
dYgdx = 1/(rhog*ug) * (omega.*w_k+mdotv*(Ydk-Yg));

dydx = [dTgdx; drhogdx; dugdx; dTddx; duddx; drddx; dYgdx];

HRR = -netProdRates(gas)'*enthalpies_RT(gas)*gasconstant*Tg;
ER = 0;
extras = [HRR, P, mdotv, ER, nd];


end
