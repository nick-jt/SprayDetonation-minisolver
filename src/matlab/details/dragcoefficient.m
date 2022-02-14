function [Cd]=dragcoefficient(y,mu,c) % state vector, viscosity, soundspeed

	% unpacking variables from state vector
	rhog 	= y(2);
	ug 	= y(3);
	ud	= y(5);
	rd	= y(6);

	%% See "Compressibility and Rarefaction Effects on Drag of a Spherical Particle" Loth

	Re = 2*Rd*(ud-ug)*rhog/mu;	% Reynolds Num
	Ma = abs(ud-ug)/c;		% Relative Mach num

	% Empirical correlation
	if (Ma<0.89)
		Gm = 1-1.525*Ma^4;
	elseif (Ma>=0.89)
		Gm = 10^(-4)*(2+8*tanh(12.77*(Ma-2.02)));
	end


	% Drag ratio
	if (Ma<1.45)
		Cm = 1/3*(5+2*tanh(3*log(Ma+0.1)));
	elseif (Ma>=1.45)
		Cm = 2.044+0.2*exp(-1.8*log(Ma/1.5)^2);
	end

	% Empirical correlation
	Hm = 1-0.258*Cm/(1+514*Gm);

	% Drag coefficient
	if (Re<0.1)
		Cd = 24/Re;
	elseif (Re<45)
		Cd = (24/Re)*(1+0.15*Re^0.687);
	elseif (Re>=45)
		Cd = Hm*(24/Re)*(1+0.15*Re^0.687)+0.42*Cm/(1+42500*Gm*Re^(-1.16));
	end

end
