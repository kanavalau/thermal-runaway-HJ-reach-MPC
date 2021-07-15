function dy_dt = ODEs_system(t,y,q,param)

Cp_r    = param(1);
Cp_c    = param(2);
rho_r   = param(3);
rho_c   = param(4);
d_H     = param(5);
k0      = param(6);
R       = param(7);
Ea      = param(8);
nA      = param(9);
nuA     = param(10);
nuB     = param(11);
V_r     = param(12);
V_c     = param(13);
A       = param(14);
Uhtc    = param(15);
% qmax    = param(16);
% dq_max  = param(17);
Tc_in   = param(18);

% Set up the equations
dy_dt = [   nuA * k0 * exp(-Ea / (R * y(3,:))) * y(1,:).^nA;                  % 1 change in conc of A
            nuB * k0 * exp(-Ea / (R * y(3,:))) * y(1,:).^nA;                  % 2 change in conc of B
            1./(rho_r*V_r*Cp_r) .* (k0 * exp(-Ea / (R * y(3,:))) * y(1,:).^nA...
                * (-d_H) * V_r - Uhtc .* A .* (y(3,:) - y(4,:)));               % 3 change in temp in reactor
            q * (Tc_in - y(4,:)) / V_c + Uhtc .* A .* (y(3,:) - y(4,:))...
            / (V_c * rho_c * Cp_c)];                                            % 4 change in temp in coolant

end