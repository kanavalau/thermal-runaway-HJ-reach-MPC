function J = objective_MPC(x,tspan_0,dt,y0_0,T_sp,param,np)

% The function simulates the system for the control and prediction horizons
% and returns the value of the objective (integral of the square of the 
% difference between the reactor temperature and the set point temperature)

% Initialise total performance function
tspan   = tspan_0;                    % Time span of 1st iteration
y0    	= y0_0;                       % Initial conditions of 1st iteration

% Perform the specified number of control steps
for n = 1:length(x)
    
    % Assign flowrate as being the controlled variable
    q      = x(1,n);
    
    % Solve the actual system
    options = odeset('Jacobian',@(t,y) Jacobian_MPC(t,y,q,T_sp,param),'RelTol',1e-5);
    sol     = ode23s(@(t,y) ODEs_MPC(t,y,q,T_sp,param), tspan, y0,options);
    
    % Extract the solution and final values
    ysol	= real(sol.y);
    
    % Update time span and initial conditions
    tspan 	= tspan + dt;
    y0    	= ysol(:,end);
    
end

% Prediction horizon
tspan       = [tspan(2), tspan(2) + (np-length(x))*dt];
y0          = ysol(:,end);

% Solve the system and extract the solution
options = odeset('Jacobian',@(t,y) Jacobian_MPC(t,y,q,T_sp,param),'RelTol',1e-5);
sol        = ode23s(@(t,y) ODEs_MPC(t,y,q,T_sp,param), tspan, y0,options);
ysol        = real(sol.y);

% Return the objective to be minimised
J        = ysol(end,end);
end