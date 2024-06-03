
%%%%% All variables are dimensional variables
clear; close all;

constants;

mkdir(expdir);

% for ns =1:length(shear_all)
for ns =1:6
    ns
    % rw_all = rw_mg(ns)
    shear = shear_all(ns)
    
    mkdir([expdir 'shear_' num2str(shear*1e3,3)]);

    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
  

    parfor i=1:length(kx_all)
        kx=kx_all(i);
        grow = NaN.*zeros(1,length(m0_all));

        m0_limit_t = kx *shear/omega + m0_limit;
        
        for j=1:length(m0_all)
            m0 = m0_all(j);    
            
            if(m0>=m0_limit_t)

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
        
                
            end

            
    
        end

     %%% Save the data
     outputname=[[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_kx' num2str(i) '.mat'];

     %%% Save outputs
     s = struct('grow',grow,'shear',shear,'kx_all',kx_all,'m0_all',m0_all, ...
         'kx',kx,'nu',nu,'kappa',kappa);  
     save(sprintf(outputname),"-fromstruct",s);

    end

end


