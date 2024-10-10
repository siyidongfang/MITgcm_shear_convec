
function [dt,Nt,tt,psi,zeta,buoy,p0,b0,z0,bq1,bq2,bq3,bq4,bq5,zq1,zq2,zq3,zq4,...
            dbdt,dzetadt,b0_wgrid,dbdz,d2bdz2,d2zetadz2,dpsidz,dUdz,U,U_wgrid] ...
            = initialize(Umax,lambda,Ptide,USEdiffusion,nu,dz,Lt,Nr)

        % %%% Estimate time step
        % CFLx = 0.8;
        % if(Umax~=0)
        %     dt_cfl = CFLx/Umax*lambda;   % The time step required to satisfy the CFL consition
        % else
        %     dt_cfl = CFLx/0.0001*lambda;
        % end
        
        dt_tide = Ptide/(72*5);       % The time step required to resolve tides
        % dt = min([dt_tide dt_cfl]);
        dt = dt_tide;
        
        % if(USEdiffusion)
        %     %%% Time step constraint based on horizontal diffusion 
        %     deltaT_Ah = 0.5*(lambda/4)^2/(4*nu);  
        %     %%% Time step constraint based on vertical diffusion
        %     deltaT_Ar = 0.5*dz^2 / (4*nu);
        %     dt = min([dt_tide dt_cfl deltaT_Ah deltaT_Ar])
        % end
        
        Nt = round(Lt/dt);
        tt = dt:dt:Nt*dt;
        
        %%%% Define variables
        psi = zeros(Nt,Nr+1);
        zeta = zeros(Nt,Nr+1);
        buoy = zeros(Nt,Nr);
        
        p0 = zeros(1,Nr+1);
        z0 = zeros(1,Nr+1);
        b0 = zeros(1,Nr);
        
        bq1 = zeros(Nt,Nr);
        bq2 = zeros(Nt,Nr);
        bq3 = zeros(Nt,Nr);
        bq4 = zeros(Nt,Nr);
        bq5 = zeros(Nt,Nr);
        dbdt = zeros(Nt,Nr);
        
        zq1 = zeros(Nt,Nr+1);
        zq2 = zeros(Nt,Nr+1);
        zq3 = zeros(Nt,Nr+1);
        zq4 = zeros(Nt,Nr+1);
        dzetadt = zeros(Nt,Nr+1);
        
        b0_wgrid = zeros(1,Nr+1);
        dbdz = zeros(1,Nr+1);
        d2bdz2 = zeros(1,Nr);
        d2zetadz2 = zeros(1,Nr+1);
        
        dpsidz = zeros(1,Nr+1);
        dUdz = zeros(1,Nr);
        U = zeros(1,Nr);
        U_wgrid = zeros(1,Nr+1);
        
        %%% Initial condition
        % buoy(1,:) = 2.0000e-23;   
        buoy(1,:) = 2.0000e-23*(randn(1,Nr)+1i*randn(1,Nr)); % Guassian white noise
        psi(1,:) = 0;
        zeta(1,:) = 0;

end
