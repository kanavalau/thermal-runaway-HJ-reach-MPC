clear;

% Add required folders to path
addpath('../','../../System_parameters/')

% Get the system parameters
system_parameters;

% Initial conditions:
cA0         = 8.;                 	% Initial concentration of reactant kmol/m^3
cB0         = 0.;                   % Initial concentration of product kmol/m^3
Tr_0        = 400;                  % Initial reactor temperature in K
Tc_0        = Tc_in;                % Initial coolant temperature in K

% Set point temperature:
T_sp     	= 410;                  % in K
     
% Specify simulation domain
ni      = 500;                      % Total number of iterations
nc      = 3;                        % Steps in control horizon
np      = 10;                       % Steps in prediction horizon
dt      = 10;                       % Control time step

% Boundedness of instantaneous change in coolant flow rate (linear inequality constraint)
A1 = zeros(2*(nc-1),nc);
A1(1:2,1)  = [-1;1];
A1(end-1:end,end)  = [1;-1];
for i = 2:nc-1
    A1(2*(i-1)-1:2*(i-1)+3-1,i) = [1;-1;-1;1];  % LHS matrix for inequality constraints
end                      
b       = dq_max*q_max*ones(2*(nc-1),1);         % RHS vector for inequality constraints

% There are no linear equality constraints
Aeq     = [];                       % Matrix coefficient for equality constraints
beq     = [];                       % Constants for equality constraints

% Control variable bounds
lb      = zeros(1,nc);              % Lower bounds on coolant flow
ub      = q_max*ones(1,nc);          % Upper bounds on coolant flow
qg      =  lb;                      % First guess of coolant flow

% Initial conditions
y0      = [cA0, cB0, Tr_0, Tc_0,0];
tspan   = [0 dt];

% Solutions vector
Sols = [];
% Control actions vector
qs = [];
% Exitflag vector
ef = zeros(1,ni);


for i = 1:ni
    
    % Calculate optimal control parameters
    [q_int, ~, exitflag, ~] = fmincon(@(x) objective_MPC(x,tspan,dt,y0,T_sp,param,np),...
            qg, A1, b, Aeq, beq, lb, ub,[]);
        
	ef(i) = exitflag;
    
    % First lower and upper bound are updated every control set to ensure
    % maximum instantaneous change is not exceeded between MPC interations
    % while linear constraints ensure that is the case between control
    % steps whithin MPC
    
    % Check if fmincon has converged:
    if exitflag == -2
        % No feasible solution - set coolant flow rate to upper bound
        q      = ub(1);
        % Reset upper and lower bounds
        if (ub(1) + dq_max * q_max < q_max)
            ub(1)   = ub(1) + dq_max * q_max;
            lb(1)   = ub(1) - 2 * dq_max * q_max;
        else
            ub(1)   = q_max;
            lb(1)   = q_max  - dq_max * q_max;
        end
    else
        
        % All other cases - set cooling to the solution found and update
        % intial guess
        q      = q_int(1,1);
        qg      = q_int(1,end-1) * ones(size(qg));
        
        % Reset lower bound
        if (q_int(1) - dq_max * q_max) < 0
            lb(1) = 0;
        else
            lb(1) = q_int(1) - dq_max * q_max;
        end
        
        % Reset upper bound
        if (q_int(1) + dq_max * q_max) > q_max
            ub(1) = q_max;
        else
            ub(1) = q_int(1) + dq_max * q_max;
        end
    end
        
    % Solve the system of differential equations
    options = odeset('Jacobian',@(t,y) Jacobian_MPC(t,y,q,T_sp,param),'RelTol',1e-5);
    sol     = ode23s(@(t,y) ODEs_MPC(t,y,q,T_sp,param), tspan, y0,options);

    % Actual solution:
    y       = real(sol.y);              % Extract the solution profile
    t       = real(sol.x);              % Extract the time
    
    % Update initial conditions and time span
    y0      = y(:,end);
    tspan   = tspan + dt;

    % Reset performance function
    y0(end) = 0;
    
    % Save the solutions obtained
    Sols    = [Sols(:,1:end-1), [t;y]];
    % Store control actions taken
    qs = [qs(1:end-1),q*ones(1,length(t))];

end


% Plot reactor temperature profile
figure(1)
plot(Sols(1,:),Sols(4,:))
xlabel('Time, $t$/s','interpreter','latex')
ylabel('Reactor temperature, $T_\mathrm{r}$/K','interpreter','latex')

% Plot coolant flow rate
figure(2)
plot(Sols(1,:),qs)
xlabel('Time, $t$/s','interpreter','latex')
ylabel('Coolant flow rate, $q$/m$^3$ s$^{-1}$','interpreter','latex')

% Plot product concentration
figure(3)
plot(Sols(1,:),Sols(3,:))
xlabel('Time, $t$/s','interpreter','latex')
ylabel('Product concentration, $c_\mathrm{B}$/kmol m$^{-3}$','interpreter','latex')

% Load the avoid set and calulate value function profile
avoid_set_data = load('../../Compute_avoid_set/avoid_set_od.mat');

alpha = zeros(length(Sols(1,:)),1);

for i = 1:length(Sols(1,:))
    alpha(i) = calc_lsf(Sols(2,i),Sols(4,i),Sols(5,i),avoid_set_data.g,avoid_set_data.data);
end

% Plot value function profile
figure(4)
plot(Sols(1,:),alpha)
xlabel('Time, $t$/s','interpreter','latex')
ylabel('Level-set function $\alpha^{*}$','interpreter','latex')

% Save data for plotting
save(['plotting_data_MPC_standard_',num2str(T_sp)],'Sols','qs','alpha')