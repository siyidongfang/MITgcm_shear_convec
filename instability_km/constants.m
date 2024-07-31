

NT1 = 20;
NT2 = 70;

pi = vpa(pi, 60);
% expdir = 'exps_new/flat_kv2e-4/';
expdir = 'exps_new/topo4_kv2e-4_pi60/';

Diffusion = true;
ConvectiveAdjustment = false;
nt_percycle = 72*30; 

N = sqrt(1)*1e-3;

topo=4;
shear_Ri0_25 = 0.0018; % 0.0017525;
shear_Ri1 = 9.7e-04;

% topo=0;
% shear_Ri0_25 = 2*N;
% shear_Ri1 = N;

Ptide = 43200;
omega = 2*pi/Ptide;
% shear_all = [0:1e-4/10:shear_Ri0_25];
% shear_all = [0:1e-4/50:shear_Ri0_25/4]; %%% for small-shear
shear_all = [0:1e-4:shear_Ri0_25]; 

Ns = length(shear_all);

h_shear = 250;
m0_rw = 2*pi/h_shear;

lam_z_all = [1:1:500];
m0_all = 2*pi./lam_z_all;

lam_x_all = [5 10:10:12000];
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




% m0_all = [m0min:0.01/10:0.6 0.61:0.01/2:1];
% kx_all = [k0min:0.001/10:0.06 0.061:0.001/2:0.1];
% lam_z_all = 2*pi./m0_all;
% lam_x_all = 2*pi./kx_all;

% lam_x_all = flip(lam_x_all);
% rw_all= 10.^([-3:0.1/2:-1 -0.95:0.01/2:0 0.1:0.1/2:2.1]);
 
% m0max = 2*pi/1;
% m0min = m0_rw;
% k0max = 2*pi/3;
% k0min = 2*pi/30000;

