
NTtide = 1;

% pi = vpa(pi, 60);
% pi_hp = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
% expdir = 'eigenvalues/eig_topo4_kv8e-5/';
expdir = 'eigenvalues/eig_flat_kv8e-5/';

Diffusion = true;
nt_percycle = 72*30; 

N = sqrt(1)*1e-3;

% topo=4;
% shear_Ri0_25 = 0.0018; % 0.0017525;
% shear_Ri1 = 9.7e-04;

topo=0;
shear_Ri0_25 = 2*N;
shear_Ri1 = N;

Ptide = 43200;
omega = 2*pi/Ptide;
% shear_all = [0:1e-4/10:shear_Ri0_25];
% shear_all = [0:1e-4/50:shear_Ri0_25/4]; %%% for small-shear
shear_all = [0:1e-4:shear_Ri0_25]; 

Ns = length(shear_all);

h_shear = 250;
m0_rw = 2*pi/h_shear;

lam_z_all = [1:1:500];
% lam_z_all = [0.01:0.01:0.09 0.1:0.1:10 10.25:0.25:50 50.5:0.5:500];
m0_all = 2*pi./lam_z_all;

lam_x_all = [5 10:10:12000];
% lam_x_all = [0.05 0.1 0.25:0.25:10 10:1:800 805:5:12000];
kx_all = 2*pi./lam_x_all;

rw_all = kx_all/m0_rw;
lam_x_real = 2*pi./(m0_rw.*rw_all);

Nk = length(kx_all);
Nm = length(m0_all);
Nrw = length(rw_all);

cs = cosd(topo);
ss = sind(topo);

kappa_const = 8e-5;
nu_const = 8e-5;

if(Diffusion)
    kappa = kappa_const;
    nu = nu_const;
else 
    kappa = 0;
    nu = 0;
end



