function status = myoutputfunc(~,y,flag,gas)

status = 0;

if strcmp(flag,'init') || strcmp(flag,'done')
else % after every timestep

    if (y==0)
	 status=1; % case returned error
    end

    set(gas,'T',y(1),'Rho',y(2),'Y',y(7:end));
    c = soundspeed(gas);
    M = y(3)/c;
    
    if M>1-1e-3
        status = 1;
    end
end
end
