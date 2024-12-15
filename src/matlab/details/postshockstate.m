function gas = postshockstate( U0, Pg0, Tg0, q, mech )

       gas = Solution(mech);
       set(gas,'T',Tg0,'P',Pg0,'X',q);
       p1 = Pg0;
       rho1 = density(gas);
       u1 = U0;
       h1 = enthalpy_mass(gas);
       
       rho2 = rho1*1.5;
       rho2_prev = 4;
       
       while abs(rho2-rho2_prev)>1e-4*rho2
            u2 = rho1*u1/rho2;
            p2 = p1+rho1*u1^2-rho2*u2^2;
            h2 = h1+0.5*u1^2-0.5*u2^2;
	    set(gas,'H',h2,'P',p2);
            rho2_prev = rho2;
            rho2 = density(gas);
       end
            
end
