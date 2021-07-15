function [c, ceq] = NLconstraints(x,tspan_0,dt,y0_0,T_sp,param,np,avoid_set_data)

data = avoid_set_data.data;
g = avoid_set_data.g;

% Initialise total performance function
tspan   = tspan_0;                    % Time span of 1st iteration
y0    	= y0_0;                       % Initial conditions of 1st iteration
    
c       = zeros(length(x)+1,1);

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
        
    c(n) = - calc_lsf(ysol(1,end),ysol(3,end),ysol(4,end),g,data);
end
    
    
tspan       = [tspan(2), tspan(2) + (np-length(x))*dt];
y0          = ysol(:,end);
    
options = odeset('Jacobian',@(t,y) Jacobian_MPC(t,y,q,T_sp,param),'RelTol',1e-5);
sol        = ode23s(@(t,y) ODEs_MPC(t,y,q,T_sp,param), tspan, y0,options);
ysol        = real(sol.y);

c(length(x)+1) = - calc_lsf(ysol(1,end),ysol(3,end),ysol(4,end),g,data);

ceq = [];
end

