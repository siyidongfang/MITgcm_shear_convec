
% clear;

Lz = 1;
dz = 0.01;
Nr = round(Lz/dz)+1;
z0 = (rand(1,Nr)-0.5)/1e30;

m1km = 1000;
lambda = 0.5*m1km;
H = 1500;
delta = H;
kx = 2*pi/(lambda/delta)
C = kx^2*dz^2-2;

p0 = zeros(1,Nr);
An = zeros(Nr-2,Nr-2); %%% Matrix
Dn = z0(2:Nr-1)'*dz^2;
An(1,1)=C;An(1,2)=1;
An(Nr-2,Nr-3)=1;An(Nr-2,Nr-2)=C;
for n=2:Nr-3
    An(n,n-1)=1;
    An(n,n)=C;
    An(n,n+1)=1;
end
det_An = det(An);
for n=1:Nr-2
    Bn = An;
    Bn(:,n) = Dn;
    det_Bn = det(Bn);
    p0(n+1) = det_Bn/det_An;
end

plot(p0);
