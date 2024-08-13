
%%%%% All variables are dimensional variables
clear; close all;

constants;
mkdir(expdir);


% for ntopo = 1:Ntopo
for ntopo = 8

    topo = topo_all(ntopo)
    cs = cosd(topo);
    ss = sind(topo);

    for nn = 1:Nn
        N = N_all(nn)

        clear shear_all;
        shear_all = [0:3*N/(Ns-1):3*N]; 


        for ns =1:Ns
            shear = shear_all(ns);
            rs = shear/omega; %%% shear over omega 
            if(omega==0)
                rs = 0;
            end
            M11 = NaN.*zeros(Nk,Nm); % zeta(T) with the initial condition (z0=1,b0=0)
            M21 = NaN.*zeros(Nk,Nm); % buoy(T) with the initial condition (z0=1,b0=0)
            M12 = NaN.*zeros(Nk,Nm); % zeta(T) with the initial condition (z0=0,b0=1)
            M22 = NaN.*zeros(Nk,Nm); % buoy(T) with the initial condition (z0=0,b0=1)
            lambda = NaN.*zeros(Nk,Nm,2); % eigenvalues of the Floquet Matrix
            grow = NaN.*zeros(Nk,Nm); % Growth rate in 1/s
        
            parfor i=1:Nk
                kx=kx_all(i);
        
                for j=1:Nm
                    m0 = m0_all(j);    
                    
                    z0 = 1;
                    b0 = 0;
                    [dt,Nt,tt,zeta,buoy,dbdt,dzetadt] = initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,b0,z0);
                    [buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t] =loop(dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta);
                    M11(i,j)=zeta(end);
                    M21(i,j)=buoy(end);
        
                    z0 = 0;
                    b0 = 1;
                    [dt,Nt,tt,zeta,buoy,dbdt,dzetadt] = initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,b0,z0);
                    [buoy,zeta,psi,www,uuu,re_buoy,re_uuu,re_www,ct,st,mz_t,angle_front,a1_t] =loop(dt,Nt,dbdt,dzetadt,omega,m0,rs,kx,shear,ss,cs,N,kappa,nu,tt,buoy,zeta);
                    M12(i,j)=zeta(end);
                    M22(i,j)=buoy(end);
        
                    % Define Floquet Matrix
                    M = [M11(i,j),M12(i,j);M21(i,j),M22(i,j)];
                    % Compute the eigenvalues and eigenvectors
                    [V, D] = eig(M);
                    % The eigenvalues are the diagonal elements of D
                    lambda(i,j,:) = diag(D);
                    % Compute the growth rate in 1/hour
                    grow(i,j) = max(real(lambda(i,j,:)));
        
                end
                
            end
        
            clear M11 M12 M21 M22 Lambda
            save([expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3,3) '_shear' num2str(shear*1e3,3) '.mat']);
        
        end
    
    end
end



