% Specific molar heat capacities:
Cp_r        = 2.33e3;               % Heat capacity of reaction mixture in J/(kg*K)
Cp_c        = 4.18e3;               % Heat capacity of cooling water J/(kg*K)

% Reaction mixture specifics:
R           = 8.314;                % Universal gas constant J/(mol*K)
Ea          = 9500*R;               % Activation energy in J/mol
d_H         = -100e6;               % Enthalpy change of formation in J/kmol
k0          = 1.5e6;                % Pre-exponential factor for reaction in 1/s
rho_r       = 950;                  % Density of reaction mixture in kg/m^3
nA          = 1.;                   % Reaction order of A
nuA         = -1.;                  % Stoichiometric coefficient for A
nuB         = 1.;                   % Stroichiometric coefficient for B

% Properties of the reactor:
V_r         = 20;                   % Volume of reactor in m^3
A           = 36;                   % Heat transfer area in m^2
V_c         = 1.4;                  % Cooling jacket volume in m^3
q_max    	= 0.030;            	% Maximum coolant flow rate in m^3/s
dq_max      = 0.1;                  % Maximum fractional step-change in coolant flow -
rho_c       = 1000;                 % Cooling water density kg/m^3
Uhtc       	= 400;                  % Heat transfer coefficient W/(m^2*K)
Tc_in       = 300;                  % Inlet coolant temperature in K

% Parameters vector
param = [Cp_r
         Cp_c
         rho_r
         rho_c
         d_H
         k0
         R
         Ea
         nA
         nuA
         nuB
         V_r
         V_c
         A
         Uhtc
         q_max
         dq_max
         Tc_in];