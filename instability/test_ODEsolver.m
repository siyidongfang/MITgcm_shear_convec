
clear;
% close all;

Lz = 1;
dz = 0.01;
Nr = round(Lz/dz)+1;
z0 = (rand(1,Nr)-0.5);
% z0 = ones(1,Nr);
% z0 = [1:Nr];
zz = 0:dz:((Nr-1)*dz);

m1km = 1000;
lambda = 0.03*m1km;
H = 500;
delta = H;
kx = 2*pi/(lambda/delta)
C = -kx^2*dz^2-2;

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

figure(1)
subplot(1,2,1)
plot(p0,zz);
subplot(1,2,2)
u=-diff(p0)/dz;
plot(u,0.5*(zz(1:end-1)+zz(2:end)))
