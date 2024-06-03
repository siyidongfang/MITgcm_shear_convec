

kappa_const = 2e-8;
nu_const = 2e-8;
NT1 = 25;
NT2 = 80;
expdir = 'parallel_flat_nu2e-8/'; %%% Flat bottom with diffusion/viscous dissipation

Diffusion = true;
ConvectiveAdjustment = false;
nt_percycle = 72*10; 

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

% m0_all = [0 m0min*[1:1/2:7.5] 0.01:0.01/10:0.6 0.61:0.01/2:1];
m0_all = [m0min:0.01/10:0.6 0.61:0.01/2:1];
kx_all = [k0min:0.001/10:0.06 0.061:0.001/2:0.1];
lam_z_all = 2*pi./m0_all;
lam_x_all = 2*pi./kx_all;

%--- constants
cs = cosd(topo);
ss = sind(topo);
    
b00 = 2.0e-23;
b0 = b00*(rand()+rand()*1i);  %%% Initial condition b(t=0)


if(Diffusion)
    kappa = kappa_const;
    nu = nu_const;
else 
    kappa = 0;
    nu = 0;
end

%--- constants
 
