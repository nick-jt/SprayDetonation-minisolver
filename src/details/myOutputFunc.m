
function status = myOutputFcn(~,y,flag,gas)

status = 0;

if flag ~= strcmp(flag,'init') | strcmp(flag,'done')
else % after every timestep
    set(gas,'T',y(1),'Rho',y(2),'Y',y(7:end));
    c = soundspeed(gas);
    M = y(3)/c;
    
    if M>1-1e-3
        status = 1;
    end
end

end