function J = objective(x,tspan_0,dt,y0_0,param,np)

% Initialise total performance function
tspan   = tspan_0;                    % Time span of 1st iteration
y0    	= y0_0;                       % Initial conditions of 1st iteration
    
    for n = 1:length(x)
        
        % Assign flowrate as being the controlled variable
        u0      = x(1,n);
        
        % Solve the actual system
        sol     = ode23s(@(t,y) ODEs_MPC(t,y,u0,param), tspan, y0);
        
        % Extract the solution and final values
        ysol	= real(sol.y);
        
        % Update time span and initial conditions
        tspan 	= tspan + dt;
        y0    	= ysol(:,end);
        
    end
    
    
    tspan       = [tspan(2), tspan(2) + (np-length(x))*dt];
    y0          = ysol(:,end);
    
    sol        = ode23s(@(t,y) ODEs_MPC(t,y,u0,param), tspan, y0);
    ysol        = real(sol.y);
    J        = ysol(end,end);
    %J        = ysol(end,end) + 10^2*norm(x(1,2:n)-x(1,1:n-1));
end