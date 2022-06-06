% Use this file to initialize the parameters for your case
function [x,y,M] = SprayDet_Dodecane
tic

addpath('details');

% Reading case details from CaseSetup.m
C = CaseSetup;
[funcflag,U0,phi,Cdw,Chw,Rd0,alpha,T0,P0,fuel,mech,satpressure,latheat,dropCv...
	,lchar,Pr,Le,Tw,rhod,Length] = C.vars{1:end};
Cvd = dropCv(0,0); % Dodecane this value doesnt matter

% Initializing case conditions
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
	fuel phi mech char(q) gas satpressure latheat dropCv};
printcase(vars,funcflag);

% Calling respective function
if (funcflag==1)
	[U0,x,y,M] = getSSvelocity(1500,2000,vars,1);
elseif (funcflag==2)
	[x,y,M] = integrator(U0,vars);
elseif (funcflag==3)
	[SS_Cdws(i),~,~,~] = bracketMethodCDW(U0,0,0.1,vars,1);
else
	error("Unknown value for funcflag in CaseSetup");
end


%%%%%%%%%Post processing
fprintf("Post-Processing Results\n");
[Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
	lam, alpha, Length, fuel, phi, mech, q, gas...
	, satpressure, latheat, dropCv] = vars{1:end};

% Initialize the initial state gas
set(gas, 'T', Tg0, 'P', Pg0, 'X', char(q));
rho0 = density(gas);

% Calculate number density (modifies gas state)
if (Rd0>0 && alpha>0)
	nd = getnumden(gas,alpha,phi,rhod,Rd0,fuel,rho0);
else
	nd = 0;
end

% Re-initialize the gas
set(gas, 'T', Tg0, 'P', Pg0, 'X', char(q));
fprintf("I.C.: T=%f P=%f rho=%f U=%f nd=%f\n",temperature(gas),pressure(gas),density(gas),U0,nd);

% Calculate Postshock Conditions
gas = postshockstate( U0, Pg0, Tg0, char(q), mech );
fprintf("Postshock: T=%f P=%f rho=%f\n",temperature(gas),pressure(gas),density(gas));
lam = thermalConductivity(gas);

Tg1 = temperature(gas);
Rhog1 = density(gas);
Ug1 = rho0*U0/Rhog1;
Td1 = Tg0;
Ud1 = U0;
Rd1 = Rd0;
Yg1 = massFractions(gas);

% Initializing variables required for postprocessing
nu0 = nd*Ud1;
vars = {Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, U0,...
	lam, alpha, Length, fuel, phi, mech, q, gas...
	, satpressure, latheat, dropCv};
extras = zeros(length(x),11);
for i = 1:length(x)
	[~,extras(i,:)] = statevectorfunction(x,y(i,:)',vars,1);
end
fprintf("here")

% Confirmed mass conservation
intmdotvdx = trapz(x,extras(:,3));
fprintf("delta(rho*u)=%f integral(mdotv,dx)=%f\n",y(end,2)*y(end,3)-rho0*U0,intmdotvdx)

if (size(x)==1)
	error("size(x)=1");
end

rmpath('details');

filename = sprintf("spray_R%.2e_T%.2f_P%.2f_Phi%.2f_%s.dat",Rd0,T0,P0,phi,fuel);
fileout = fopen(filename,'w');
fprintf(fileout,"X[m] Tg[K] Pg[Pa] Rhog[kg/m^3] Rd[m] Yf Ug[m/s] HRR ER DiffTimescale[s] VaporTimescale[s] MinTimescale[s] CO_timescale[s] OH_timescale[s] Fuel_timescale[s] Temp_timescale[s]\n");
for i=1:length(x)
	fprintf(fileout,"%.10f %f %f %f %e %e %f %e %f %e %f %e %e %e %e %e\n",...
		x(i),y(i,1),extras(i,2),y(i,2),y(i,6),y(i,6+29),y(i,3),... %7
		extras(i,1),extras(i,4),extras(i,5),extras(i,6),extras(i,7),... %5
		extras(i,8),extras(i,9),extras(i,10),extras(i,11)); %4
end
toc

end
