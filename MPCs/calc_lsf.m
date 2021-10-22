function c = calc_lsf(cA,Tr,Tc,g,data)

% The function takes system state and value function data and computes the
% value function for the state using spline interpolation

Y = linspace(g.min(1),g.max(1),g.N(1));
X = linspace(g.min(2),g.max(2),g.N(2));
Z = linspace(g.min(3),g.max(3),g.N(3));

if Tr<g.min(2)
    Tr = g.min(2);
elseif Tr>g.max(2)
    Tr = g.max(2);
end

if Tc<g.min(3)
    Tc = g.min(3);
elseif Tc>g.max(3)
    Tc = g.max(3);
end

if cA<g.min(1)
    cA = g.min(1);
elseif cA>g.max(1)
    cA = g.max(1);
end

c = interp3(X,Y,Z,data,Tr,cA,Tc,'spline');

end

