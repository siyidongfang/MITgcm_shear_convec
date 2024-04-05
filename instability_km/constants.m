
b00 = 1e-70;
b0 = b00*(rand()+rand()*1i);  %%% Initial condition b(t=0)
Lt = NTtide*Ptide; 
dt = Ptide/nt_percycle;


Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;

kappa = 1e-7;
nu = 1e-6;

psi = zeros(1,Nt);
zeta = zeros(1,Nt);
buoy = zeros(1,Nt);
dbdt = zeros(1,Nt);
dzetadt = zeros(1,Nt);

%%% Initial condition
buoy(1) = b0;
psi(1) = 0;
zeta(1) = 0;

cs = cosd(topo);
ss = sind(topo);
