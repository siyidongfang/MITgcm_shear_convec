clear

phi1=0;
omega = 2*pi/43200;
omega1 = 8*omega;
omega2 = 1*omega;
a1 = 1;

syms  phi2 phi3 a2 a3 

eq1 = a1*cos(omega1*pi+phi1)-a2*cos(omega2*pi+phi2)==0;
eq2 = omega1*a1*sin(omega1*pi+phi1)-omega2*a2*sin(omega2*pi+phi2)==0;

eq3 = a2*cos(omega2*2*pi+phi2)-a3*cos(omega1*2*pi+phi3)==0;
eq4 = omega2*a2*sin(omega2*2*pi+phi2)-omega1*a3*sin(omega1*2*pi+phi3)==0;

S = solve(eq1,eq2,eq3,eq4);

n=4;
phi2 = double(S.phi2(n));
phi3 = double(S.phi3(n));
a2 =double(S.a2(n));
a3 =double(S.a3(n));

% a2*cos(omega2*2*pi+phi2)-1
% 
% abs(a3)-a1

