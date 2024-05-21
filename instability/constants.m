
fontsize = 16;

kx = 2*pi/lambda;

t1hour = 3600;
m1km = 1000;

N = 1e-3;
topo = 0;
Ptide = 43200;
omega = 2*pi/Ptide;
NTtide = 50;
Lt = NTtide*Ptide; 

Hmax = 250;
Umax = Hmax * Shear;

dz = 1;      
Nr = round(Hmax/dz);
zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
zz_wgrid = 0:dz:((Nr)*dz);

if(USEdiffusion)
    % nu = 2e-6; 
    % kappa = 2e-6;
    nu = 2e-4; 
    kappa = 2e-4;
else
    nu = 0;
    kappa = 0;
end
Pr = nu/kappa;

CFLx = 0.5;
if(Umax~=0)
    dt_cfl = CFLx/Umax*lambda;   % The time step required to satisfy the CFL consition
else
    dt_cfl = CFLx/0.0001*lambda;
end

dt_tide = Ptide/(72*10)/12;       % The time step required to resolve tides
dt = min([dt_tide dt_cfl]);

if(USEdiffusion)
    %%% Time step constraint based on horizontal diffusion 
    deltaT_Ah = 0.5*(lambda/4)^2/(4*nu)   
    %%% Time step constraint based on vertical diffusion
    deltaT_Ar = 0.5*dz^2 / (4*nu)
    dt = min([dt_tide dt_cfl deltaT_Ah deltaT_Ar])
end

Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;

%%%% Define variables
psi = zeros(Nt,Nr+1);
zeta = zeros(Nt,Nr+1);
buoy = zeros(Nt,Nr);

p0 = zeros(1,Nr+1);
z0 = zeros(1,Nr+1);
b0 = zeros(1,Nr);

bq1 = zeros(Nt,Nr);
bq2 = zeros(Nt,Nr);
bq3 = zeros(Nt,Nr);
bq4 = zeros(Nt,Nr);
bq5 = zeros(Nt,Nr);
dbdt = zeros(Nt,Nr);

zq1 = zeros(Nt,Nr+1);
zq2 = zeros(Nt,Nr+1);
zq3 = zeros(Nt,Nr+1);
zq4 = zeros(Nt,Nr+1);
dzetadt = zeros(Nt,Nr+1);

b0_wgrid = zeros(1,Nr+1);
dbdz = zeros(1,Nr+1);
d2bdz2 = zeros(1,Nr);
d2zetadz2 = zeros(1,Nr+1);

dpsidz = zeros(1,Nr+1);
dUdz = zeros(1,Nr);
U = zeros(1,Nr);
U_wgrid = zeros(1,Nr+1);

%%% Initial condition
buoy(1,:) = 2.0000e-23;   
psi(1,:) = 0;
zeta(1,:) = 0;
