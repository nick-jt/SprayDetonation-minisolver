
function [SS_Cdw,t,y,M] = bracketMethodCDW(V,CDWlow,CDWhigh,vars,fileout)

iter = 0;
while abs(CDWlow-CDWhigh)>1e-4
    
    iter = iter + 1;
    CDW = (CDWlow+CDWhigh)/2;
    fprintf(fileout,"Iter %d CDW=%f",iter,CDW);
    vars{3} = CDW;
    
    [t,y,M] = integrator(V,vars);
    
    
    if (t(end)<vars{16} || (t(end)>=vars{16} && M(end)>0.99))
        CDWlow = CDW;
    else
        CDWhigh = CDW;
    end
    
    fprintf(fileout,"(t_end=%f,M_end=%f)\n",t(end),M(end));
    
end

SS_Cdw = (CDWlow+CDWhigh)/2;
end
