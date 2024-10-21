
%%%%% All variables are dimensional variables
clear; close all;
NT2 = 70;

N = 1e-3;

Diffusion = false;
ConvectiveAdjustment = false;
nt_percycle = 72*30; 

Ptide = 43200;
omega = 2*pi/Ptide;

% topo=0;
% load('../figures/fig4/Ri_flat.mat')
% shear_Ri0_25 = 2*N;
% shear_Ri1 = N;
% Nri = 693; %% flat, 1/Ri=3
% % rj = 2.5000e-102; %% flat, 1/Ri=3
% rj  =  1e-150;
% % Nri = 401; %% flat, 1/Ri=1
% % rj =14.56; %% flat, 1/Ri=1

topo=4;
load('../figures/fig4/Ri_topo4_new.mat')
shear_Ri0_25 = 0.0018; % 0.0017525;
shear_Ri1 = 9.7e-04;
% Nri = 650; %% topo4, 1/Ri_min=2.5
% rj = cotd(4)-shear_calc_Ri(2923)/omega; %%%m0/k0

% Nri = 700 %% topo4, 1/Ri_min = 3
% rj = cotd(4)-shear_calc_Ri(3148)/omega %%%m0/k0

Nri = 702 %% topo4, 1/Ri_min = 3
rj = cotd(4)-shear_calc_Ri(3160)/omega %%%m0/k0

% % Nri = 699; %% topo4, 1/Ri_min=3
% % rj = 3.3220; %% topo4, 1/Ri_min=3
% % Nri = 432; %% topo4, 1/Ri_min=1
% % rj = 6.2820; %% topo4, 1/Ri_min=1


max_shear = shear_Ri0_25/2*3;
Ns = 1200;
shear_all = [0:max_shear/(Ns-1):max_shear]; 


h_shear = 250;
m0_rw = 2*pi/h_shear;

lam_z_all = [1:1:500];
m0_all = 2*pi./lam_z_all;

lam_x_all = [1e-100 5 10:10:12000];
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



for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end



inverseRi = 1./Ri_km(Nri)


%%

m0 = m0_rw;

grow_rw = zeros(Ns,Nrw);

for ns =Nri
    ns
    shear = shear_all(ns)
    
    rs = shear/omega; %%% shear over omega 

    for j=1

        % rw = rw_all(j);
        rw = 1./rj;
        kx = m0*rw;

        grow = NaN.*zeros(1,Nrw);
                    
     
        NTtide = NT2*2;

        [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
        initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

        [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        =loop(grow,j,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);

    end
end

save('budget/topo4_output_Ri3_new.mat');

