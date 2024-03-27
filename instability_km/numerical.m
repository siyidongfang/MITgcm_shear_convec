%%%%% All variables are dimensional variables

clear all;

%%% Define variables
N = 1e-3;
shear = 1.3e-3;
Ptide = 43200;
omega = 2*pi/Ptide;
kx = 1;
mz = 1;
rw = kx/mz; %%% ratio of the wavenumbers kx/mz
rs = shear/omega; %%% shear over omega
b0 = 1e-5;  %%% Initial condition b(t=0)
topo =2;

NTtide = 5;
dt = 60;
Lt = NTtide*43200; 
Nt = Lt/dt;
tt = dt:dt:Nt*dt;

psi = zeros(1,Nt);
zeta = zeros(1,Nt);
buoy = zeros(1,Nt);
dbdt = zeros(1,Nt);
dzetadt = zeros(1,Nt);

cs = cosd(topo);
ss = sind(topo);

%%% Initial condition
buoy(1) = b0;
psi(1) = 0;
zeta(1) = 0;

%%% Start the loop
for o=1:Nt-1
    %%% Fourth-order Runge-Kutta method %%%
    t0 = tt(o);
    b0 = buoy(o);
    z0 = zeta(o);
    tendency;
    k_1b = dbdt(o);
    k_1z = dzetadt(o);
    % Euler forward predictor advancing dt/2:
    b_2 = buoy(o)+0.5*dt*k_1b;
    z_2 = zeta(o)+0.5*dt*k_1z;
    t0 = tt(o)+dt/2;
    b0 = b_2;
    z0 = z_2;
    tendency;
    k_2b = dbdt(o);
    k_2z = dzetadt(o);
    % Euler backward corrector advancing dt/2:
    b_3 = buoy(o)+0.5*dt*k_2b;
    z_3 = zeta(o)+0.5*dt*k_2z;
    t0 = tt(o)+dt/2;
    b0 = b_3;
    z0 = z_3;
    tendency;
    k_3b = dbdt(o);
    k_3z = dzetadt(o);
    % Mid-point predictor advancing dt:
    b_4 = buoy(o)+dt*k_3b;
    z_4 = zeta(o)+dt*k_3z;
    t0 = tt(o)+dt;
    b0 = b_4;
    z0 = z_4;
    tendency;
    k_4b = dbdt(o);
    k_4z = dzetadt(o);

    % Simpson rule corrector advancing dt:
    buoy(o+1) = buoy(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
    zeta(o+1) = zeta(o) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;

end

ct = cos(omega*tt);
st = sin(omega*tt);
a1t = -(kx^2+mz^2+2*kx*mz*rs*st+kx^2*rs^2*st.^2);
a2t = 1i*kx*(2*rs*ss*st-cs)*N^2 - 1i*mz*N^2*ss;
a3t = 1i*kx*(cs-rs*ss*st) - 1i*mz*ss;

psi = zeta./a1t;
www = 1i*kx*psi;
uuu = -1i*mz*psi-1i*kx*psi*rs.*st;

pe = buoy.^2;
ke = real(0.5*(uuu.^2+www.^2));
kew = real(0.5*(www.^2));

figure(1)
clf
plot(tt/43200,buoy)
hold on;
plot(tt/43200,real(www)/400)

figure(2)
clf;
plot(tt/43200,pe);
hold on;
plot(tt/43200,kew/1e5);
ylim([0 max(pe)/10])





