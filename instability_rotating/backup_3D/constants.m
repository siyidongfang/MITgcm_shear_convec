

NT1 = 10;
NT2 = 40;

expdir = 'exps_test/';
N = 1e-3;
f0 = 1.2e-4;

Diffusion = false;
ConvectiveAdjustment = false;
nt_percycle = 72*30; 

topo=4;
shear_Ri0_25 = 0.0018; % 0.0017525;
shear_Ri1 = 9.7e-04;

% topo=0;
% shear_Ri0_25 = 2*N;
% shear_Ri1 = N;

Ptide = 43200;
omega = 2*pi/Ptide;
max_shear = shear_Ri0_25/2*3;
% Ns = 1200;
Ns = 14;
shear_all = [0:max_shear/(Ns-1):max_shear]; 

% shear_all = [0:1e-4/50:shear_Ri0_25/4]; %%% for small-shear
% shear_all = [0:1e-4:shear_Ri0_25]; 
% Ns = length(shear_all);

h_shear = 250;
m0_rw = 2*pi/h_shear;

Lz_all = [-500:2:500];
m0_all = 2*pi./Lz_all;

Lx_all = [-12000:25:12000];
k0_all = 2*pi./Lx_all;

Ly_all = [-12000:25:12000];
l0_all = 2*pi./Ly_all;

rw_all = k0_all/m0_rw;
lam_x_real = 2*pi./(m0_rw.*rw_all);

Nk = length(k0_all);
Nl = length(l0_all);
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


