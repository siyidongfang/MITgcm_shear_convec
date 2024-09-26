
NT1 = 80;
% NT2 = 70;


expdir = 'exps_rotating/noV_flat_N1e-3';
N = 1e-3;
f = 1.2e-4;

Diffusion = false;
ConvectiveAdjustment = false;
nt_percycle = 72*30; 


% topo=4;
% shear_Ri0_25 = 0.0018; % 0.0017525;
% shear_Ri1 = 9.7e-04;

topo=0;
shear_Ri0_25 = 2*N;
shear_Ri1 = N;

Ptide = 43200;
omega = 2*pi/Ptide;
max_shear = shear_Ri0_25/2*2.5;
Ns = 1000;
shear_all = [0:max_shear/(Ns-1):max_shear]; 

% shear_all = [0:1e-4/50:shear_Ri0_25/4]; %%% for small-shear
% shear_all = [0:1e-4:shear_Ri0_25]; 
% Ns = length(shear_all);

h_shear = 250;
m0_rw = 2*pi/h_shear;

lam_z_all = [1:1:500];
m0_all = 2*pi./lam_z_all;

% lam_x_all = [-6000:10:6000];
lam_x_all = [-6000:5:6000];
kx_all = 2*pi./lam_x_all;

rw_all = kx_all/m0_rw;
lam_x_real = 2*pi./(m0_rw.*rw_all);

Nk = length(kx_all);
Nm = length(m0_all);
Nrw = length(rw_all);

cs = cosd(topo);
ss = sind(topo);
    
b00 = 2.0e-70;
b0 = b00*(rand()+rand()*1i);  %%% Initial condition b(t=0)

kappa_const = 2e-4;
nu_const = 2e-4;

if(Diffusion)
    kappa = kappa_const;
    nu = nu_const;
else 
    kappa = 0;
    nu = 0;
end

