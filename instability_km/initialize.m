
function [dt,Nt,tt,psi,zeta,buoy,dbdt,dzetadt,dbdz_vert,dBdz_vert,dB0dz_vert,dbtotaldz_vert] = initialize(shear,h_shear,kx,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,ConvectiveAdjustment,b0)


        lambda = 2*pi/kx;
        dz = 2*pi/m0;

        % Umax = shear*h_shear;
        % %%% Estimate time step
        % CFLx = 0.5;
        % if(Umax~=0)
        %     dt_cfl = CFLx/Umax*lambda;   % The time step required to satisfy the CFL consition
        % else
        %     dt_cfl = CFLx/0.0001*lambda;
        % end

        dt_tide = Ptide/nt_percycle;       % The time step required to resolve tides
        % dt = min([dt_tide dt_cfl]);
        dt = dt_tide;

        if(Diffusion)
            %%% Time step constraint based on horizontal diffusion
            deltaT_Ah = 0.5*(lambda/4)^2/(4*nu);
            %%% Time step constraint based on vertical diffusion
            deltaT_Ar = 0.5*dz^2 / (4*nu);
            % dt = min([dt_tide dt_cfl deltaT_Ah deltaT_Ar])
            dt = min([dt_tide deltaT_Ah deltaT_Ar]);
        end


    Lt = NTtide*Ptide;
    dt = Ptide/nt_percycle;
    
    Nt = round(Lt/dt);
    tt = dt:dt:Nt*dt;
    
    if(omega==0)
        Nt = 1e3;
        dt = NTtide*Ptide/Nt;
        tt = dt:dt:Nt*dt;
    end

    psi = zeros(1,Nt);
    zeta = zeros(1,Nt);
    buoy = zeros(1,Nt);
    dbdt = zeros(1,Nt);
    dzetadt = zeros(1,Nt);
    
    if(ConvectiveAdjustment)
        dbdz_vert = zeros(1,Nt);
        dBdz_vert = zeros(1,Nt);
        dB0dz_vert = zeros(1,Nt);
        dbtotaldz_vert = zeros(1,Nt);
    else
        dbdz_vert =[];
        dBdz_vert =[];
        dB0dz_vert=[];
        dbtotaldz_vert=[];
    end
    
    %%% Initial condition
    buoy(1) = b0;
    psi(1) = 0;
    zeta(1) = 0;
    

end
