classdef CaseSetup %% Here you can define parameters for your case
	properties(Constant)

	% Flag defines what you would like your case to do:
	% 1 = iterate the eigenvalue solution to obtain SS detonation velocity
	% 2 = take provided velocity and calculate detonation structure (may not be SS ZND solution)
	% 3 = take provided velocity and obtain steady state drag coefficient

	funcflag = 1;

	% Most common parameters
	U0 	= 1790; 	% Detonation velocity (treated differently based on funcflag)
	phi     = 1.0;		% equivalence ratio
	Cdw     = 0.00;		% wall drag coefficient
	Chw     = 0.00;		% wall heat loss coefficient
	Rd0     = 5e-6; 	% drop Radius
	alpha   = 1;		% droplet loading (kg liquid fuel/kg total fuel)
	T0      = 420;		% initial temp (K)
	P0      = 101325;	% initial pressure (pa)
	fuel    = 'NC12H26';     % Change below parameters if this changes
	mech    = 'Dodecane.cti'; % Cantera chemical mechanism

	% Droplet related parameters
	satpressure = nan;      % Will be calculated in code (Dodecane only)
	latheat = @(Td,wf) 260628; % J/kg
	dropCv = @(Td,w) 2766.0;

	% Parameters
	lchar   = 3.81*0.01/4;
	Pr      = 1;
	Le      = 1;
	Tw      = CaseSetup.T0;
	rhod    = 750;
	Length  = 0.2;

	vars = {CaseSetup.funcflag,CaseSetup.U0,CaseSetup.phi,CaseSetup.Cdw,CaseSetup.Chw,CaseSetup.Rd0,CaseSetup.alpha,CaseSetup.T0...
		,CaseSetup.P0,CaseSetup.fuel,CaseSetup.mech,CaseSetup.satpressure,CaseSetup.latheat,CaseSetup.dropCv...
		,CaseSetup.lchar,CaseSetup.Pr,CaseSetup.Le,CaseSetup.Tw,CaseSetup.rhod,CaseSetup.Length};

	end
end
