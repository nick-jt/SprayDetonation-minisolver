function nd = getnumden(gas,alpha,phi,rhod,rd,fuel,rho0)

C_count = nAtoms(gas,fuel,'C');
H_count = nAtoms(gas,fuel,'H');
a       = C_count + 0.25 * H_count; % Air stoichiometric coefficient
q       = fuel + ":" + string(phi) ...
            + ", O2:" + string(a) ...
            + ", N2:" + string(a*3.76); % TODO
set(gas,'X',char(q));
Y = massFractions(gas);
fuel_index = speciesIndex(gas,fuel);
Ynf = Y; Ynf(fuel_index) = [];

MF = alpha*Y(fuel_index) / (Y(fuel_index)*(1-alpha) + sum(Ynf));
nd = rho0*MF/(rhod*4/3*pi*rd^3);

end
