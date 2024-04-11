
Lt = NTtide*Ptide; 
dt = Ptide/nt_percycle;

Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;


if(omega==0)
    Nt = 1e3;
    dt = NTtide*Ptide/Nt;
    tt = dt:dt:Nt*dt;
end

b00 = 1e-5;
b0 = b00*(rand()+rand()*1i);  %%% Initial condition b(t=0)

kappa_const = 1e-7;
nu_const = 1e-6;

% kappa_const = 1e-6;
% nu_const = 1e-5;

if(Diffusion)
    kappa = kappa_const;
    nu = nu_const;
else 
    kappa = 0;
    nu = 0;
end

psi = zeros(1,Nt);
zeta = zeros(1,Nt);
buoy = zeros(1,Nt);
dbdt = zeros(1,Nt);
dzetadt = zeros(1,Nt);

if(ConvectiveAdjustment)
dbdz_vert = zeros(1,Nt);
dBdz_vert = zeros(1,Nt);
dB0dz_vert = zeros(1,Nt);
dbtotaldz_vert = zeros(1,Nt);
end

%%% Initial condition
buoy(1) = b0;
psi(1) = 0;
zeta(1) = 0;

cs = cosd(topo);
ss = sind(topo);
