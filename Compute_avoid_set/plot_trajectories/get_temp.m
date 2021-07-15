function Tr = get_temp(cA,Tc)

load('../avoid_set_od.mat','data','g');

index_cA = (cA-g.min(1))*(g.N(1)-1)/(g.max(1)-g.min(1))+1;

index_Tc = (Tc-g.min(3))*(g.N(3)-1)/(g.max(3)-g.min(3))+1;
        
x0 = floor(index_cA);
x1 = ceil(index_cA);


z0 = floor(index_Tc);
z1 = ceil(index_Tc);

xd = index_cA - x0;

zd = index_Tc - z0;

c0 = data(x0,:,:)*(1-xd) + data(x1,:,:)*xd;
c1 = c0(1,:,z0)*(1-zd) + c0(1,:,z1)*zd;

i_min = length(nonzeros(c1>0));

i_av = i_min + c1(i_min)/(c1(i_min)+abs(c1(i_min+1)));

Tr = (i_av-1)*(g.max(2)-g.min(2))/(g.N(2)-1) + g.min(2);
end