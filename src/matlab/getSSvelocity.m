function [SS_Cdw,t,y,M] = getSSvelocity(Vlow,Vhigh,vars,fileout)

iter = 0;
while abs(Vlow-Vhigh)>1e-4

    iter = iter + 1;
    Vmid = (Vlow+Vhigh)/2;
    fprintf(fileout,"Iter %d V=%f",iter,Vmid);

    [t,y,M]=integrator(Vmid,vars);


    if (t(end)<vars{16} || (t(end)>=vars{16} && M(end)>0.99))
        Vlow = Vmid;
    else
        Vhigh = Vmid;
    end

    fprintf(fileout,"(t_end=%f,M_end=%f)\n",t(end),M(end));

end

SS_Cdw = (Vlow+Vhigh)/2;
end
