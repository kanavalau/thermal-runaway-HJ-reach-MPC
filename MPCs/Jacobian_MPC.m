function df_dy = Jacobian_MPC(t,y,q,T_sp,param)

% The function computes the Jacobian corresponding to the ODEs_MPC.m
% differential equations

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
% Tc_in   = param(18);

df_dy = zeros(length(y(:,1)),length(y(:,1)));

df_dy(1,1) = nuA * k0 * exp(-Ea ./ (R * y(3,:)));
df_dy(1,2) = 0;
df_dy(1,3) = nuA * k0 * exp(-Ea ./ (R * y(3,:))) .* y(1,:).^nA .* (Ea ./ (R * y(3,:).^2));
df_dy(1,4) = 0;
df_dy(1,5) = 0;

df_dy(2,1) = nuB * k0 * exp(-Ea ./ (R * y(3,:)));
df_dy(2,2) = 0;
df_dy(2,3) = nuB * k0 * exp(-Ea ./ (R * y(3,:))) .* y(1,:).^nA .* (Ea ./ (R * y(3,:).^2));
df_dy(2,4) = 0;
df_dy(2,5) = 0;

df_dy(3,1) = 1./(rho_r*Cp_r) .* (k0 * exp(-Ea ./ (R * y(3,:))) * (-d_H));
df_dy(3,2) = 0;
df_dy(3,3) = 1./(rho_r*V_r*Cp_r) .* (k0 * exp(-Ea ./ (R * y(3,:))) .* y(1,:).^nA .* (Ea ./ (R * y(3,:).^2))...
                * (-d_H) * V_r - Uhtc .* A);
df_dy(3,4) = 1./(rho_r*V_r*Cp_r) .* (Uhtc .* A);
df_dy(3,5) = 0;

df_dy(4,1) = 0;
df_dy(4,2) = 0;
df_dy(4,3) = 1./(rho_c*V_c*Cp_c) .* (Uhtc .* A);
df_dy(4,4) = - q/V_c - 1./(rho_c*V_c*Cp_c) .* (Uhtc .* A);
df_dy(4,5) = 0;

df_dy(5,1) = 0;
df_dy(5,2) = 0;
df_dy(5,3) = 2 * (y(3,:)-T_sp);
df_dy(5,4) = 0;
df_dy(5,5) = 0;

end

