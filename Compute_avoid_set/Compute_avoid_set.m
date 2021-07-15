function [ data, g, data0 ] = Compute_avoid_set(accuracy)
%
% This code uses "A Toolbox of Level Set Methods." Ian M. Mitchell.
% UBC Department of Computer Science Technical Report TR-2007-11 (June 2007).
% https://www.cs.ubc.ca/~mitchell/ToolboxLS/
%
% Based on ToolboxLS/Examples/Reachability/air3D.m
%
% Input Parameters:
%
%   accuracy: Controls the order of approximations.
%     'low': Use odeCFL1 and upwindFirstFirst.
%     'medium': Use odeCFL2 and upwindFirstENO2 (default).
%     'high': Use odeCFL3 and upwindFirstENO3.
%     'veryHigh': Use odeCFL3 and upwindFirstWENO5.
%
% Output Parameters:
%
%   data: Implicit surface function at t_max.
%
%   g: Grid structure on which data was computed.
%
%   data0: Implicit surface function at t_0.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
% run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 5000;                    % End time.
plotSteps = 5;                  % How many intermediate plots to produce?
t0 = 0;                         % Start time.
singleStep = 0;                 % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'local';

%---------------------------------------------------------------------------
% Problem Parameters
% Add required folders to path
addpath('../System_parameters/')

% Get the system parameters
system_parameters;

% Expected maximum reactant concentration:
cA0         = 8.;                 	% kmol/m^3

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the 3D reachable set.
displayType = 'surface';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Visualize the angular dimension a little bigger.
%aspectRatio = [ 1 1 1 ];

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

%---------------------------------------------------------------------------
% Approximately how many grid cells?
Nx = 50;

% Create the grid.
% Define the domain
g.dim = 3;
g.min = [  0; 400; Tc_in];
g.max = [ cA0; 450; Tc_in+50];
g.bdry = { @addGhostExtrapolate; @addGhostExtrapolate; @addGhostExtrapolate };

% Discretise the grid
g.N = [ Nx; 4*Nx; Nx/2 ];
g = processGrid(g);

%---------------------------------------------------------------------------
% Create initial conditions (hyperplane at Tr=450 with normal in the negative direction).
data = shapeHyperplane(g,[0;-1;0],[4;450;325]);
data0 = data;

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @HamFunc;
schemeData.partialFunc = @PartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.Tc_in = Tc_in;
schemeData.Cp_r = Cp_r;
schemeData.Cp_c = Cp_c;
schemeData.R = R;
schemeData.Ea = Ea;
schemeData.d_H = d_H;
schemeData.k0 = k0;
schemeData.rho_r = rho_r;
schemeData.nA = nA;
schemeData.nuA = nuA;
schemeData.nuB = nuB;
schemeData.V_r = V_r;
schemeData.A = A;
schemeData.V_c = V_c;
schemeData.q_max = q_max;
schemeData.dq_max = dq_max;
schemeData.rho_c = rho_c;
schemeData.Uhtc = Uhtc;



%---------------------------------------------------------------------------
% Choose degree of dissipation.

switch(dissType)
 case 'global'
  schemeData.dissFunc = @artificialDissipationGLF;
 case 'local'
  schemeData.dissFunc = @artificialDissipationLLF;
 case 'locallocal'
  schemeData.dissFunc = @artificialDissipationLLLF;
 otherwise
  error('Unknown dissipation function %s', dissFunc);
end

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.
%   The Lax-Friedrichs approximation scheme MUST already be completely set up.
innerFunc = schemeFunc;
innerData = schemeData;
clear schemeFunc schemeData;

% Wrap the true Hamiltonian inside the term approximation restriction routine.
schemeFunc = @termRestrictUpdate;
schemeData.innerFunc = innerFunc;
schemeData.innerData = innerData;
schemeData.positive = 0;

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);

camlight right;  camlight left;
hold on;
axis(g.axis);
%daspect(aspectRatio);
drawnow;

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t, y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Get correct figure, and remember its current view.
  figure(f);
  [ view_az, view_el ] = view;

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);

  % Restore view.
  view(view_az, view_el);
  
end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);

save('avoid_set_od','data','data0','g')
savefig('avoid_set_plot')

%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = HamFunc(t, data, deriv, schemeData)

% Function defining the Hamiltonian

checkStructureFields(schemeData, 'grid', 'k0', 'Ea', ...
                                 'R', 'd_H','V_r','Uhtc','A','rho_r',...
                                 'Cp_r','V_c','rho_c','Cp_c','q_max',...
                                 'Tc_in');

grid = schemeData.grid;

hamValue = -(-deriv{1}.*schemeData.k0.*exp(-schemeData.Ea./(schemeData.R*grid.xs{2})).*grid.xs{1} ...
                +deriv{2}.*(schemeData.k0*exp(-schemeData.Ea./(schemeData.R*grid.xs{2})).*grid.xs{1}*(-schemeData.d_H)*schemeData.V_r-schemeData.Uhtc*schemeData.A*(grid.xs{2}-grid.xs{3}))/(schemeData.V_r*schemeData.rho_r*schemeData.Cp_r)...
                +deriv{3}.*schemeData.Uhtc*schemeData.A.*(grid.xs{2}-grid.xs{3})/(schemeData.V_c*schemeData.rho_c*schemeData.Cp_c)...
                +schemeData.q_max*deriv{3}.*(schemeData.Tc_in-grid.xs{3})/(schemeData.V_c).*(deriv{3}.*(schemeData.Tc_in-grid.xs{3})/(schemeData.V_c)>0));


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = PartialFunc(t, data, derivMin, derivMax, schemeData, dim)

% Function defining the partial differentials of the Hamiltonian w.r.t.
% derivatives of the value function
% Required by the Lax-Friedrichs method

checkStructureFields(schemeData, 'grid', 'k0', 'Ea', ...
                                 'R', 'd_H','V_r','Uhtc','A','rho_r',...
                                 'Cp_r','V_c','rho_c','Cp_c','q_max',...
                                 'Tc_in');

grid = schemeData.grid;

switch dim
  case 1
    alpha = abs(-schemeData.k0*exp(-schemeData.Ea./(schemeData.R*grid.xs{2})).*grid.xs{1});

  case 2
    alpha = abs((schemeData.k0*exp(-schemeData.Ea./(schemeData.R*grid.xs{2})).*grid.xs{1}*(-schemeData.d_H)*schemeData.V_r-schemeData.Uhtc*schemeData.A*(grid.xs{2}-grid.xs{3}))/(schemeData.V_r*schemeData.rho_r*schemeData.Cp_r));

  case 3
    alpha = abs(schemeData.Uhtc*schemeData.A*(grid.xs{2}-grid.xs{3})/(schemeData.V_c*schemeData.rho_c*schemeData.Cp_c))+abs((schemeData.Tc_in-grid.xs{3})/(schemeData.V_c)*schemeData.q_max);
    
  otherwise
    error([ 'Partials for the exothermic batch reaction' ...
            ' only exist in dimensions 1-3' ]);
end