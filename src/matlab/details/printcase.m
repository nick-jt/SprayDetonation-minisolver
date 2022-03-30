function printcase(vars,funcflag)
[Tg0, Pg0, Cdw, Chw, Rd0, lchar, Pr, Le, Tw, Cvd, rhod, nu0, D,...
    lam, alpha, Length, fuel, phi, mech, q, gas] = vars{1:end};

fprintf('SprayDetonation-minisolver\n');
fprintf('	Copyright (c) 2021 Nicolas Tricard\n');
fprintf('	Please cite:\n'); 
fprintf('		"One dimensional steady-state modeling\n');
fprintf('	of spray detonations considering loss effects"\n');
fprintf('	Nicolas Tricard, Xinyu Zhao\n');
fprintf('		"Heterogeneous effects on the propagation\n');
fprintf('	and quenching of spray detonations" Tianfeng Lu\n');
fprintf('	and Chung Law\n');
fprintf('\nCase Setup:\n');
fprintf('	Fuel		  = %s\n',fuel);
fprintf('	Det Velocity      = %.2f [m/s]\n',D);
fprintf('	Droplet radius    = %.3e [m]\n',Rd0);
fprintf('	Equivalence Ratio = %.2f [-]\n',phi);
fprintf('	Drag coefficient  = %.2f [-]\n',Cdw);
fprintf('	Heat loss coeff   = %.2f [-]\n',Chw);
fprintf('	T0		  = %.2f [K]\n',Tg0);
fprintf('	P0		  = %.2e [Pa]\n',Pg0);

if (funcflag==1)
	fprintf("Determining steady state velocity\n");
elseif (funcflag==2)
	fprintf("Integrating with provided velocity %f\n",D);
elseif (funcflag==3)
	fprintf("Calculating SS drag coefficient\n");
else
	error("Unknown value for funcflag in CaseSetup");
end

end
