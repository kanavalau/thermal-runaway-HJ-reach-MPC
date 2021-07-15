clear;

% Add required folders to path
addpath('../../System_parameters/')

% Get the system parameters
system_parameters;

% Initial conditions:
cA0         = 8.;                 	% Initial concentration of reactant kmol/m^3
cB0         = 0.;                   % Initial concentration of product kmol/m^3
Tc_0        = Tc_in;                % Initial coolant temperature in K

% Open the avoid set plot
open('../avoid_set_plot.fig')
% Make it transparent
h = get(gca,'Children');
set(h(1),'FaceAlpha',0.5)
hold on
     
% Get temperature on the level set surface
Tr_0 = get_temp(cA0,Tc_0);
     
% Initial conditions outside the avoid set
y0      = [cA0, cB0, Tr_0-1, Tc_0];
tspan   = [0 5000];

% Solve the system of differential equations
options = odeset('Jacobian',@(t,y) Jacobian_system(t,y,q_max,param),'RelTol',1e-5);
sol     = ode23s(@(t,y) ODEs_system(t,y,q_max,param), tspan, y0, options);

% solution:
y       = real(sol.y);              % Extract the solution profile
t       = real(sol.x);              % Extract the time

% Plot the trajectory (outside the avoid set)
plot3(y(1,:),y(3,:),y(4,:),'b','LineWidth',2)
% Plot the projection below the trajectory
c300 = 300*ones(length(y(4,:)),1);
plot3(y(1,:),y(3,:),c300,'Color',[0.2 0.2 0.2])
% Specify the number of vertical lines
N_vl = 6;
% Plot the vertical lines
plot3([y(1,1),y(1,1)],[y(3,1),y(3,1)],[y(4,1),300],'Color',[0.2 0.2 0.2])
for i=1:N_vl-1
    [~,j] = min(abs(y(3,:)-(400+(max(y(3,:))-400)/(N_vl-1)*i)));
    plot3([y(1,j),y(1,j)],[y(3,j),y(3,j)],[y(4,j),300],'Color',[0.2 0.2 0.2])
end

% Initial conditions inside the avoid set
y0      = [cA0, cB0, Tr_0+1, Tc_0];
tspan   = [0 5000];

% Solve the system of differential equations
options = odeset('Jacobian',@(t,y) Jacobian_system(t,y,q_max,param),'RelTol',1e-5);
sol     = ode15s(@(t,y) ODEs_system(t,y,q_max,param), tspan, y0, options);

% solution:
y       = real(sol.y);              % Extract the solution profile
t       = real(sol.x);              % Extract the time

% Plot the trajectory (insider the avoid set)
plot3(y(1,:),y(3,:),y(4,:),'r','LineWidth',2)
% Plot the projection below the trajectory
c300 = 300*ones(length(y(4,:)),1);
plot3(y(1,:),y(3,:),c300,'Color',[0.2 0.2 0.2])
% Specify the number of vertical lines
N_vl = 6;
% Plot the vertical lines
plot3([y(1,1),y(1,1)],[y(3,1),y(3,1)],[y(4,1),300],'Color',[0.2 0.2 0.2])
for i=1:N_vl-1
    [~,j] = min(abs(y(3,:)-(min(y(3,:))+(450-min(y(3,:)))/(N_vl-1)*i)));
    plot3([y(1,j),y(1,j)],[y(3,j),y(3,j)],[y(4,j),300],'Color',[0.2 0.2 0.2])
end

% Add axis labels
zlabel('Coolant temperature, $T_\mathrm{c}$/K','Interpreter','latex');
ylabel('Reactor temperature, $T_\mathrm{r}$/K','Interpreter','latex');
xlabel('Concentration of A, $c_\mathrm{A}$/kmol m$^{-3}$',...
    'Interpreter','latex');

% Set the view angle
view(gca,[71.4 17.9824817518248]);
% Save figure
savefig('avoid_set_with_traj')