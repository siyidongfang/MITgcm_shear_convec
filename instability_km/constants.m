

NT1 = 40;
NT2 = 80;
expdir = 'parallel_flat_rw/'; %%% Flat bottom with diffusion/viscous dissipation

Diffusion = false;
ConvectiveAdjustment = false;
nt_percycle = 72*50; 

topo=0;
N = sqrt(1)*1e-3;
Ptide = 43200;
omega = 2*pi/Ptide;
shear_Ri0_25 = 2*N;
shear_Ri1 = N;
shear_all = [0:1e-4:shear_Ri0_25];

h_shear = 1000;
m0_limit = 2*pi/h_shear;

m0max = 2*pi/1;
m0min = m0_limit;
k0max = 2*pi/3;
k0min = 2*pi/30000;

m0_all = [m0min:0.01/10:0.6 0.61:0.01/2:1];
kx_all = [k0min:0.001/10:0.06 0.061:0.001/2:0.1];
lam_z_all = 2*pi./m0_all;
lam_x_all = 2*pi./kx_all;

Nk = length(kx_all);
Nm = length(m0_all);

m0_rw = 2*pi/500;
rw_all= 10.^([-3:0.1/2:-1 -0.95:0.01/2:0 0.1:0.1/2:2.1]);
lam_x_real = 2*pi./(m0_rw.*rw_all);

Nrw = length(rw_all);

cs = cosd(topo);
ss = sind(topo);
    
% b00 = 2.0e-23;
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

 
