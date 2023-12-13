
clear;close all;

kx = 2;

Nr = 1000;

zspan = [0 1];
% av = [1:100];
av = -rand(1,100)+1i*rand(1,100);  
% av = [1 2 3 4 10 20 40 60 100 200 300 500 2 6000 100000 40 20 10 3 2 1 0]
zeta = @(z) interp1(linspace(min(zspan),max(zspan),numel(av)), av, z);    

%%% Form Initial Guess
xmesh = linspace(0,1,Nr/10);
solinit = bvpinit(xmesh, [0.1 0]);
%%% Solve the BVP using the bvp4c solver.
sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta), @bcfun, solinit);

%%% Solve the BVP a second time using a different initial guess for the solution.
solinit = bvpinit(xmesh, [3 0]);
sol2 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta), @bcfun, solinit);


%%% Compare Solutions
plot(sol1.x,real(sol1.y(1,:)),'-o',sol2.x,real(sol2.y(1,:)),'-o')
title('BVP with Different Solutions That Depend on the Initial Guess')
xlabel('x')
ylabel('y')
legend('Solution 1','Solution 2')

dpsidz = real(sol1.y(2,:));
psi = real(sol1.y(1,:));

xxx  = sol1.x;
aaa = diff(psi)./diff(xxx);
figure()
plot((xxx(1:end-1)+xxx(2:end))/2,aaa)
hold on;
plot(xxx,dpsidz)

figure()
plot(real(av))

%%% Code Equation
function dydz = bvpfun(z,y,kx,zeta)
    dydz = [y(2); kx^2*y(1)+zeta(z)];
end

%%% Code Boundary Conditions
function res = bcfun(ya,yb)
    res = [ya(1)
           yb(1)];
end

