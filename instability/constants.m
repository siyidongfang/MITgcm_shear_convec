
    fontsize = 16;
    
    t1hour = 3600;
    m1km = 1000;
    
    N = 1e-3;
    topo = 0;
    Ptide = 43200;
    omega = 2*pi/Ptide;
    NTtide = 13/100;
    Lt = NTtide*Ptide; 
    
    h_shear = 250;
    dz = 2;      
    
    if(useLinearShear)
        Hmax = h_shear;
        Nr = round(Hmax/dz);
        zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
        zz_wgrid = 0:dz:((Nr)*dz);
    end
    
    if(useTanhShear)
        Hmax = h_shear+250;
        Nr = round(Hmax/dz);
        zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
        zz_wgrid = 0:dz:((Nr)*dz);
    end
    

    if(USEdiffusion)
        % nu = 2e-6; 
        % kappa = 2e-6;
        nu = 2e-4; 
        kappa = 2e-4;
    else
        nu = 0;
        kappa = 0;
    end
    Pr = nu/kappa;
    
    CFLx = 0.5;
    if(Umax~=0)
        dt_cfl = CFLx/Umax*lambda;   % The time step required to satisfy the CFL consition
    else
        dt_cfl = CFLx/0.0001*lambda;
    end
    
    dt_tide = Ptide/(72*2);       % The time step required to resolve tides
    dt = min([dt_tide dt_cfl]);
    
    if(USEdiffusion)
        %%% Time step constraint based on horizontal diffusion 
        deltaT_Ah = 0.5*(lambda/4)^2/(4*nu)   
        %%% Time step constraint based on vertical diffusion
        deltaT_Ar = 0.5*dz^2 / (4*nu)
        dt = min([dt_tide dt_cfl deltaT_Ah deltaT_Ar])
    end
    
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
    buoy(1,:) = 2.0000e-23;   
    psi(1,:) = 0;
    zeta(1,:) = 0;



    %%% Colors
    black = [0 0 0];
    verydarkgray = [0.25 0.25 0.25];
    darkgray = [0.5 0.5 0.5];
    gray = [0.7 0.7 0.7];
    boxcolor = [0.85 0.85 0.85];
    lightgray = [249 249 249]/255;
    
    
    red = [0.6350 0.0780 0.1840];
    lightred = [249 102 102]/255;
    orange = [0.8500 0.3250 0.0980];
    coral = [255 127 80]/255;
    pink = [255 153 204]/255;
    
    yellow = [0.9290 0.6940 0.1250];
    gold = [255 215 0]/255;
    brown = [153 102 51]/255;
    
    blue = [0 0.4470 0.7410];
    lightblue = [0.3010 0.7450 0.9330];
    purple = [0.4940 0.1840 0.5560];
    
    green = [0.4660 0.6740 0.1880];
    green2 = [0 153 0]/255;
    seagreen = [46 139 87]/255;
    olive = [107 142 35]/255;
    darkgreen = [21 71 52]/255;
    
    cyan = [0 255 255]/255;
    
    lightpurple = [204 153 255]/255;
    purple = [153 51 255]/255;
    darkpurple = [102 0 204]/255;
    brown1 = [153 76 0]/255;
    brown2 = [255 178 103]/255;


