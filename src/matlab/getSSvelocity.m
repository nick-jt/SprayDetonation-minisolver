function [SS_velocity,t,y,M] = getSSvelocity(Vlow,Vhigh,vars,fileout)

iter = 0;
L = inf;
while abs(Vlow-Vhigh)>1e-2 || L<vars{16}

    iter = iter + 1;
    Vmid = (Vlow+Vhigh)/2;
    fprintf(fileout,"Iter %d V=%f",iter,Vmid);

    [t,y,M]=integrator(Vmid,vars);


    if (t(end)<vars{16} || (t(end)>=vars{16} && M(end)>0.999))
        Vlow = Vmid;
    else
        Vhigh = Vmid;
    end

    L = t(end);

    fprintf(fileout,"	(t_end=%f,M_end=%f)\n",t(end),M(end));

end

SS_velocity = (Vlow+Vhigh)/2;
end
