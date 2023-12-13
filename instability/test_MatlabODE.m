
clear;close all;

kx = 2e-3;

zspan = [0 1];
av = -rand(1,100)/0.8+1i*rand(1,100);                                                         
zeta = @(z) interp1(linspace(min(zspan),max(zspan),numel(av)), av, z);    

[Z,Y] = ode45(@(z,y)diffvar(z,y,kx,zeta), zspan, [0;0]);

figure
plot(Z, real(Y))
grid


dpsidz = Y(:,2);
psi = Y(:,1);

aaa = diff(psi)./diff(Z);
figure()
plot((Z(1:end-1)+Z(2:end))/2,aaa)

figure()
plot(real(av))

function dydz = diffvar(z,y,kx,zeta)
    dydz = [y(2); kx^2*y(1)+zeta(z)];
end




% pv0 = [real(psi_init); imag(psi_init)];
% zspan = [0 1];
% [kx,zz,yv] = ode45(@imaginaryODE, zspan, pv0);
% y = yv(:,1) + 1i*yv(:,2);
% 
% 
% function fv = imaginaryODE(kx,zz,yv)
%   % Construct y from the real and imaginary components
%   y = yv(1) + 1i*yv(2);            
% 
%   % Evaluate the function
%   d2psidz2 = complexf(kx,psi,zeta);             
% 
%   % Return real and imaginary in separate components
%   fv = [real(d2psidz2); imag(d2psidz2)]; 
% end   
% 
% function d2psidz2 = complexf(kx,psi,zeta)
% 
%   d2psidz2 = kx^2*psi+zeta;
% 
% end