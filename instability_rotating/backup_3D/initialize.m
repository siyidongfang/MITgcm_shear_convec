
function [dt,Nt,tt,uu,vv,ww,bb,pp,dudt,dvdt,dwdt,dbdt] = initialize(k0,m0,Diffusion,nu,NTtide,Ptide,nt_percycle,omega,b0)


        lambda = 2*pi/k0;
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

    uu = zeros(1,Nt);
    vv = zeros(1,Nt);
    ww = zeros(1,Nt);
    bb = zeros(1,Nt);
    pp = zeros(1,Nt);
    dudt = zeros(1,Nt);
    dvdt = zeros(1,Nt);
    dwdt = zeros(1,Nt);
    dbdt = zeros(1,Nt);
    
    
    %%% Initial condition
    bb(1) = b0;
    uu(1) = 0;
    vv(1) = 0;
    ww(1) = 0;
    pp(1) = 0;
    

end
