
%%%%% All variables are dimensional variables
clear; close all;

constants;
mkdir(expdir);

% for ns =1:length(shear_all)
for ns =15
    ns
    shear = shear_all(ns)
    
    rs = shear/omega; %%% shear over omega 
    if(omega==0)
        rs = 0;
    end
  
    grow_s = zeros(length(k0_all),length(l0_all),length(m0_all));

    for i=1:length(kx_all)
        i
        k0=k0_all(i);
        grow = NaN.*zeros(length(l0_all),length(m0_all));
        
        parfor j=1:length(l0_all)
            l0=l0_all(i);

            for k=1:length(m0_all)
                m0 = m0_all(k);    
                
                NTtide = NT1;
                if(omega==0)
                    NTtide = 1/rw/shear/Ptide*10;
                end
                [dt,Nt,tt,uu,vv,ww,bb,pp,dudt,dvdt,dwdt,dbdt] = initialize(k0,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,b0);

                [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
                =loop(grow,k,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,k0,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
        
                if(grow(j,k)>0)
                    NTtide = NT2;
        
                    [dt,Nt,tt,uu,vv,ww,bb,pp,dudt,dvdt,dwdt,dbdt] = initialize(k0,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,b0);
                    
                    [grow,buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot,pp,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert]...
                    =loop(grow,k,NTtide,kappa_const,dt,Nt,dbdt,dzetadt,omega,m0,rs,k0,shear,ss,cs,N,kappa,nu,tt,buoy,zeta,Diffusion,ConvectiveAdjustment,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert);
                end
    
    
            end
        end
        grow_s(i,:)=grow;

    end

    save([expdir 'shear' num2str(shear*1e3,3) '_output.mat']);

end


