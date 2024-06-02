
%%%%% All variables are dimensional variables
clear; close all;
Diffusion = true;
ConvectiveAdjustment = false;
nt_percycle = 72*10; 


%%%%%% exps_flat_diff %%%%%%
expdir = 'parallel_flat_nu2e-4/'; %%% Flat bottom with diffusion/viscous dissipation
topo=0;
N = sqrt(1)*1e-3;
Ptide = 43200;
omega = 2*pi/Ptide;
shear_Ri0_25 = 2*N;
shear_Ri1 = N;
shear_all = [0:1e-4:shear_Ri0_25];
m0max = 2*pi/1;
m0min = 2*pi/5000;
k0max = 2*pi/3;
k0min = 2*pi/100000;
m0_all = [0 m0min*[1:1/2:7.5] 0.01:0.01/10:0.6 0.61:0.01/2:1];
kx_all = [0 k0min*[1:1/2:15] 0.001:0.001/10:0.06 0.061:0.001/2:0.1];
lam_z_all = 2*pi./m0_all;
lam_x_all = 2*pi./kx_all;


%--- constants
cs = cosd(topo);
ss = sind(topo);
    
b00 = 2.0e-23;
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
%--- constants
 

mkdir(expdir);

% for ns =1:length(shear_all)
for ns =2:3
    ns
    % rw_all = rw_mg(ns)
    shear = shear_all(ns)
    
    mkdir([expdir 'shear_' num2str(shear*1e3,3)]);

    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
   
    parfor m=1:length(m0_all)
        m
	    m0 = m0_all(m);
        grow = NaN.*zeros(1,length(kx_all));

    for i=1:length(kx_all)
        kx=kx_all(i);
        
        NTtide = 30;
        if(omega==0)
            NTtide = 1/rw/shear/Ptide*10;
        end
        [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

        [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        =loop(grow,i,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);

        if(grow(i)>0)
            NTtide = 100;

            [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

            [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
            =loop(grow,i,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        end

        if(grow(i)>0 && grow(i)<0.1)
            NTtide = 400;
            
            [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

            [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
            =loop(grow,i,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        end

    end

     %%% Save the data
     outputname=[[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m) '.mat'];

     %%% Save outputs
     s = struct('grow',grow,'shear',shear,'kx_all',kx_all,'m0_all',m0_all, ...
         'm0',m0,'dt',dt,'nu',nu,'kappa',kappa);  
     save(sprintf(outputname),"-fromstruct",s);

    end

    
end


