
    addpath ../analysis/colormaps/    
    Shear_parm = ([0.1:0.1:2.8])*1e-3;    %%% flat
    % Shear_parm = ([0.1:0.1:2.0 2.07])*1e-3; %%% topo4
    lambda_parm = [50:25:700 750:50:1000 1200:200:2000 2400:400:3200 4000:1000:8000 10000:2000:12000]; 
    % lambda_parm = [round(10.^[1.7:0.05:3 3.1:0.1:3.4 3.6 3.8 4]/10)*10];
    lambda_parm = flip(lambda_parm);
    lambda_parm = [lambda_parm round(10.^[1.6:-0.1:0.7]) 3];

    % lambda_parm = lambda_parm(42:82);

    % Ptide_parm = [0.5:0.5:5 10000]*43200;
    % topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    % N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4 5 6 7 8 9 10]*1e-3;
    
    useLinearShear = true;
    useTanhShear = false;
    USEdiffusion = true;  %%% Add diffusion/dissipation
 
    fontsize = 16;
    
    t1hour = 3600;
    m1km = 1000;
    
    N = 1e-3;
    topo = 0;
    Ptide = 43200;
    omega = 2*pi/Ptide;
    NTtide = 10;
    Lt = NTtide*Ptide; 
    
    h_shear = 250;
        
    
    if(useLinearShear)
        dz = 1;  
        Hmax = h_shear;
        Nr = round(Hmax/dz);
        zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
        zz_wgrid = 0:dz:((Nr)*dz);
    end
    
    if(useTanhShear)
        dz = 2; 
        Hmax = h_shear+250;
        Nr = round(Hmax/dz);
        zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
        zz_wgrid = 0:dz:((Nr)*dz);
    end

    zspan = [0 Hmax];

    if(USEdiffusion)
        nu = 1e-6; 
        kappa = 1e-6;
        % nu = 2e-4; 
        % kappa = 2e-4;
    else
        nu = 0;
        kappa = 0;
    end
    Pr = nu/kappa;
    

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


