function c = calc_lsf(cA,Tr,Tc,g,data)

% The function takes system state and value function data and computes the
% value function for the state using linear interpolation

index(1) = (cA-g.min(1))*(g.N(1)-1)/(g.max(1)-g.min(1))+1;

if index(1)<1
    index(1)=1;
end

if Tr > g.max(2)
	index(2) = g.N(2);
elseif Tr < g.min(2)
	index(2) = 1;
else
    index(2) = (Tr-g.min(2))*(g.N(2)-1)/(g.max(2)-g.min(2))+1;
end

if Tc > g.max(3)
	index(3) = g.N(3);
elseif Tr < g.min(3)
	index(3) = 1;
else
    index(3) = (Tc-g.min(3))*(g.N(3)-1)/(g.max(3)-g.min(3))+1;
end
        
x0 = floor(index(1));
x1 = ceil(index(1));
y0 = floor(index(2));
y1 = ceil(index(2));
z0 = floor(index(3));
z1 = ceil(index(3));

xd = index(1) - x0;
yd = index(2) - y0;
zd = index(3) - z0;

c00 = data(x0,y0,z0)*(1-xd) + data(x1,y0,z0)*xd;
c01 = data(x0,y0,z1)*(1-xd) + data(x1,y0,z1)*xd;
c10 = data(x0,y1,z0)*(1-xd) + data(x1,y1,z0)*xd;
c11 = data(x0,y1,z1)*(1-xd) + data(x1,y1,z1)*xd;

c0 = c00*(1-yd) + c10*yd;
c1 = c01*(1-yd) + c11*yd;

c = (c0*(1-zd) + c1*zd);

end
