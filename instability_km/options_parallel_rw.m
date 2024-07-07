
%%%%% All variables are dimensional variables
clear; close all;

constants;

mkdir(expdir);

m0 = m0_rw;

for ns =1:length(shear_all)
    ns
    shear = shear_all(ns)
    
    mkdir([expdir 'shear_' num2str(shear*1e3,3)]);

    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
  
    parfor j=1:Nrw
        % for j=1
        % load("grow_rw.mat")
        % rw = rw_max(ns);

        rw = rw_all(j);
        kx = m0*rw;

        grow = NaN.*zeros(1,Nrw);
                    
        NTtide = NT1;
        if(omega==0)
            NTtide = 1/rw/shear/Ptide*10;
        end
        [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

        [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        =loop(grow,j,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);

        if(grow(j)>0)
            NTtide = NT2;

            [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

            [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
            =loop(grow,j,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        end

        if(grow(j)>0 && grow(j)<0.1)
            NTtide = NT2*2;

            [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

            [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
            =loop(grow,j,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        end

         %%% Save the data
         outputname=[[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(j) '.mat'];

         %%% Save outputs
         s = struct('grow',grow,'shear',shear,'rw_all',rw_all,'m0',m0);  
         save(sprintf(outputname),"-fromstruct",s);
            
    end

end

save([expdir '.mat']);
