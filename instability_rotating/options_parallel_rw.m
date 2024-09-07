
%%%%% All variables are dimensional variables
clear; close all;

constants

m0 = m0_rw;

grow_rw = zeros(Ns,Nrw);

for ns =1:Ns
    ns
    shear = shear_all(ns)
    
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
        [dt,Nt,tt,psi,zeta,buoy,vvel,dvdt,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
            initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);

        [grow,vvel,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        =loop(grow,j,NTtide,kappa_const,dt,Nt,dvdt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,f,tt,buoy,vvel,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);

        % if(grow_ke(j)>0)
        %     NTtide = NT1;
        % 
        %     [dt,Nt,tt,psi,zeta,buoy,vvel,dvdt,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
        %     initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);
        % 
        %     [grow,vvel,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        %     =loop(grow,j,NTtide,kappa_const,dt,Nt,dvdt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,f,tt,buoy,vvel,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        % end

        % if(grow(j)>0 && grow(j)<0.15)
        %     NTtide = NT2*2;
        % 
        %     [dt,Nt,tt,psi,zeta,buoy,vvel,dvdt,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = ...
        %     initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0);
        % 
        %     [grow,vvel,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
        %     =loop(grow,j,NTtide,kappa_const,dt,Nt,dvdt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,f,tt,buoy,vvel,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        % end

         grow_rw(ns,j)=grow(j);
    end
end

save([expdir 'output.mat']);

