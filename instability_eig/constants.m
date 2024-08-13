
NTtide = 1;

expdir = 'shear_methieu_noDiff/';

Diffusion = false;
nt_percycle = 72*30; 

% Ptide_all = [1 2]*43200;
% topo_all = [0:2:20];
% N_all = [0 0.1 0.5 1:2:11]*1e-3;

topo_all = 4;
N_all = 1e-3;

Ptide = 43200;
omega = 2*pi/Ptide;

Nn = length(N_all);
Ntopo = length(topo_all);

h_shear = 250;
m0_rw = 2*pi/h_shear;

Ns = 251; % length of shear_all

kappa_const = 1.4e-7;
nu_const = 1e-6;

if(Diffusion)
    kappa = kappa_const;
    nu = nu_const;

    lam_z_all = [-500:1:500];
    lam_x_all = [-12000:10:12000];
else 
    kappa = 0;
    nu = 0;

    lam_z_all = h_shear;
    lam_x_all = [-8000:10:8000];
end

m0_all = 2*pi./lam_z_all;
kx_all = 2*pi./lam_x_all;
Nk = length(kx_all);
Nm = length(m0_all);



% lam_z_all = [0.01:0.01:0.09 0.1:0.1:10 10.25:0.25:50 50.5:0.5:500];
% lam_x_all = [0.05 0.1 0.25:0.25:10 10:1:800 805:5:12000];
% pi = vpa(pi, 60);
% pi_hp = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;



