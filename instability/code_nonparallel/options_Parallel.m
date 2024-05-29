
clear all;
close all;
addpath ../analysis/colormaps/

FigureIsVisible = 'off';

topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4 5 6 7 8 9 10]*1e-3;
Shear_parm = ([0:0.1:2.0])*1e-3;
% lambda_parm = [5 10 50:50:350 400:50:1000 1200:200:5000 6000:1000:20000 30000:10000:100000];
lambda_parm = [round(10.^[1.7:0.05:3 3.1:0.1:3.4 3.6 3.8 4]/10)*10];
lambda_parm = flip(lambda_parm);
lambda_parm = [lambda_parm round(10.^[1.6:-0.1:0.5])];
Ptide_parm = [0.5:0.5:5 10000]*43200;

exppath = 'exps_linear_dz0.5/';


parfor Nexp_lambda =1:length(lambda_parm)

        lambda = lambda_parm(Nexp_lambda);
        expfolder = [exppath 'lambda' num2str(lambda) '/']
        mkdir(expfolder); 

      for Nexp_shear =1:length(Shear_parm)
        Shear = Shear_parm(Nexp_shear)
        USEdiffusion = true;  %%% Add diffusion/dissipation
        

        %-------% Begin of the script Constants.m
                fontsize = 16;
                
                kx = 2*pi/lambda;
                t1hour = 3600;
                m1km = 1000;
                
                N = 1e-3;
                topo = 0;
                Ptide = 43200;
                omega = 2*pi/Ptide;
                NTtide = 10;
                Lt = NTtide*Ptide; 
                
                useLinearShear = true;
                useTanhShear = false;
                
                h_shear = 250;
                dz = 0.5;      
                
                if(useLinearShear)
                    Hmax = h_shear;
                    Nr = round(Hmax/dz);
                    zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
                    zz_wgrid = 0:dz:((Nr)*dz);
                    Atide = Shear*zz;
                    Atide_wgrid = Shear*zz_wgrid;
                    Umax = h_shear * Shear;
                end
                
                if(useTanhShear)
                    Hmax = h_shear+250;
                    Nr = round(Hmax/dz);
                    zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography   
                    zz_wgrid = 0:dz:((Nr)*dz);
                    
                    %%% Zero velocity at sea floor
                    Atide = h_shear*Shear *(1+ tanh( (zz  -Hmax/2) / (h_shear/2) )) /2;
                    Atide_wgrid = h_shear*Shear *(1+ tanh( (zz_wgrid  -Hmax/2) / (h_shear/2) )) /2;
                    Umax = h_shear * Shear;
                
                    % %%% Zero velocity at center
                    % Atide = h_shear*Shear *(tanh( (zz  -Hmax/2) / (h_shear/2) )) /2;
                    % Atide_wgrid = h_shear*Shear *(tanh( (zz_wgrid  -Hmax/2) / (h_shear/2) )) /2;
                    % Umax = h_shear * Shear /2;
                
                end
                
                dAdz = diff(Atide_wgrid)/dz;
                
                
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
                
                dt_tide = Ptide/(72*5);       % The time step required to resolve tides
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
        %-------% End of the script Constants.m

        expdir = [expfolder 'topo' num2str(topo) '_H' num2str(Hmax) ...
            '_N' num2str(N) '_S' num2str(Shear) ...
            '_lambda' num2str(lambda) '/'];
        outputname = [expdir 'output.mat'];
        mkdir(expdir);

        useRK4 = true;        %%% Use Tourth-order Runge-Kutta method
        useAB3 = false;       %%% Use Third-order Adams-Bashforth method

        noBQ2 = false;
        noBQ3 = false;
        noBQ4 = false;
        noZQ2 = false;
        noZQ3 = false;

        hydrostatic = false;
        % numerical;




%%% Background tidal velocity
Utide =cos(tt*omega)'.*Atide;
% Utide = repmat(cos(tt*omega)',[1 length(Atide)])...
%     .*repmat(Atide,[length(tt) 1])/U0;

% h=figure(1);
% set(h,'Visible', FigureIsVisible);clf;
% plot(Atide,zz,'LineWidth',2);
% title('Atide (m/s)');set(gca,'Fontsize',fontsize)
% grid on;grid minor;
% % saveas(h,[expdir 'fig1.png'])
% 
% h=figure(2);
% set(h,'Visible', FigureIsVisible);clf;
% plot(diff(Atide)./diff(zz),0.5*(zz(1:end-1)+zz(2:end)),'LineWidth',2);
% title('Shear (1/s)');set(gca,'Fontsize',fontsize)
% grid on;grid minor;
% % saveas(h,[expdir 'fig2.png'])
% 
% h=figure(3);
% set(h,'Visible', FigureIsVisible);clf;
% pcolor(tt/3600,zz,Utide');shading flat;colormap redblue; colorbar;
% title('Atide (m/s)');xlabel('Time (hours)');set(gca,'Fontsize',fontsize)
% saveas(h,[expdir 'fig3.png'])

% close all;

%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Fourth-order Runge-Kutta or Third-order Adams-Bashforth
%%%% or Euler forward scheme for time advancement

zspan = [0 Hmax];

%%%%%%%%%%%% B.C.-1 %%%%%%%%%%%%
zeta(1,1) = 0; zeta(1,Nr+1) = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for o=1:Nt-1

    t0 = tt(o);
    b0 = buoy(o,:);
    z0 = zeta(o,:);
    
    % loop;


    %%%%%%%%%%%% B.C.-6 %%%%%%%%%%%%
    z0(1) = 0; z0(Nr+1) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Solve the boundary value problem using the bvp4c solver.
    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%% Form Initial Guess
        % options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NMax', 20000);
        % p0_guess = real(p0);
        % solinit.y(1,:) = p0_guess;
        % solinit.y(2,2:end-1) = (p0_guess(3:end)-p0_guess(1:end-2))/dz/2; %%% Centered difference
        % solinit.y(2,1)=(p0_guess(2)-p0_guess(1))/dz;
        % solinit.y(2,end)=(p0_guess(end)-p0_guess(end-1))/dz;
    xmesh = zz_wgrid;
    solinit = bvpinit(xmesh, [0 0]);
    options = bvpset('NMax', 10000);
    if(~hydrostatic)
        %%% non-hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit,options);
    else
        %%% hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun_hydro(z,y,zeta0), @bcfun, solinit,options);
    end
    psi0 = sol1.y(1,:);
    zz0 = sol1.x;
    if(length(psi0)>length(p0))
        p0 = interp1(zz0,psi0,zz_wgrid);
        % warning('length(psi0)>length(p0)')
    else
        p0 = psi0;
    end


    %%%%%%%%%%%% B.C.-7 %%%%%%%%%%%%
    %%% Impermeable (w=0) At the ocean bottom
    p0(1) = 0;
    %%% Impermeable (w=0) At the upper boundary 
    p0(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    U = cos(omega*t0)*Atide;
    U_wgrid = cos(omega*t0)*Atide_wgrid;

    % dUdz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz;

    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz; %%% on w-grid

    %%%%%%%%%%%% B.C.-8 %%%%%%%%%%%%
    % %%% No buoyancy flux, Ocean bottom, when diffusivity is not zero
    % %%% dbdz+dB/dz+dB0/dz = 0 ==>  dbdz = -(dB/dz+dB0/dz)
    % dbdz(1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    % %%% No buoyancy flux, Upper boundary, when diffusivity is not zero
    % dbdz(Nr+1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    dbdz(1) = 0;
    dbdz(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    %%%%%%%%%%%% B.C.-9 %%%%%%%%%%%%
    b0_wgrid(1) = 2*b0(1)-b0_wgrid(2);
    b0_wgrid(Nr+1) = 2*b0(Nr)-b0_wgrid(Nr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d2bdz2 = (dbdz(2:Nr+1)-dbdz(1:Nr))/dz; %%% 2nd-order centered difference

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*U.*b0;
    bq2(o,:) = -1i*kx*p0_ugrid*N^2*cosd(topo);
    bq3(o,:) = dpsidz*N^2*sind(topo);
    bq4(o,:) = 1i*kx*p0_ugrid.*dAdz/omega*N^2*sind(topo)*sin(omega*t0);
    bq5(o,:) = kappa*(d2bdz2-kx^2.*b0);

    if(noBQ2)
        bq2(o,:) = 0;
    end
    if(noBQ3)
        bq3(o,:) = 0;
    end
    if(noBQ4)
        bq4(o,:) = 0;
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);

    % for m = 2:Nr
    %     d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    % end
    d2zetadz2(2:Nr) = (z0(1:Nr-1)-2*z0(2:Nr)+z0(3:Nr+1))/dz^2;

    %%%%%%%%%%%% B.C.-11 %%%%%%%%%%%%
    %%% Linear extrapolation
    d2zetadz2 = interp1(zz_wgrid(2:Nr),d2zetadz2(2:Nr),zz_wgrid,'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zq1(o,:) = -1i*kx*U_wgrid.*z0;
    zq2(o,:) = 1i*kx*b0_wgrid*cosd(topo);
    zq3(o,:) = -dbdz*sind(topo);
    zq4(o,:) = nu*(d2zetadz2-kx^2.*z0);

    if(noZQ2)
        zq2(o,:) = 0;
    end
    if(noZQ3)
        zq3(o,:) = 0;
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);


    psi(o,:) = p0;

    if(useRK4)
        

    %%% Fourth-order Runge-Kutta method %%%
        k_1b = dbdt(o,:);
        k_1z = dzetadt(o,:);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o,:)+0.5*dt*k_1b;
        z_2 = zeta(o,:)+0.5*dt*k_1z;
    %%%%%%%%%%%% B.C.-3 %%%%%%%%%%%%
        % z_2(1) = delta/Hshear*cos(t0); %%% No total stress
        z_2(1) = 0; z_2(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
            

    %%%%%%%%%%%% B.C.-6 %%%%%%%%%%%%
    z0(1) = 0; z0(Nr+1) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Solve the boundary value problem using the bvp4c solver.
    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%% Form Initial Guess
        % options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NMax', 20000);
        % p0_guess = real(p0);
        % solinit.y(1,:) = p0_guess;
        % solinit.y(2,2:end-1) = (p0_guess(3:end)-p0_guess(1:end-2))/dz/2; %%% Centered difference
        % solinit.y(2,1)=(p0_guess(2)-p0_guess(1))/dz;
        % solinit.y(2,end)=(p0_guess(end)-p0_guess(end-1))/dz;
    xmesh = zz_wgrid;
    solinit = bvpinit(xmesh, [0 0]);
    options = bvpset('NMax', 10000);
    if(~hydrostatic)
        %%% non-hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit,options);
    else
        %%% hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun_hydro(z,y,zeta0), @bcfun, solinit,options);
    end
    psi0 = sol1.y(1,:);
    zz0 = sol1.x;
    if(length(psi0)>length(p0))
        p0 = interp1(zz0,psi0,zz_wgrid);
        % warning('length(psi0)>length(p0)')
    else
        p0 = psi0;
    end


    %%%%%%%%%%%% B.C.-7 %%%%%%%%%%%%
    %%% Impermeable (w=0) At the ocean bottom
    p0(1) = 0;
    %%% Impermeable (w=0) At the upper boundary 
    p0(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    U = cos(omega*t0)*Atide;
    U_wgrid = cos(omega*t0)*Atide_wgrid;

    % dUdz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz;

    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz; %%% on w-grid

    %%%%%%%%%%%% B.C.-8 %%%%%%%%%%%%
    % %%% No buoyancy flux, Ocean bottom, when diffusivity is not zero
    % %%% dbdz+dB/dz+dB0/dz = 0 ==>  dbdz = -(dB/dz+dB0/dz)
    % dbdz(1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    % %%% No buoyancy flux, Upper boundary, when diffusivity is not zero
    % dbdz(Nr+1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    dbdz(1) = 0;
    dbdz(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    %%%%%%%%%%%% B.C.-9 %%%%%%%%%%%%
    b0_wgrid(1) = 2*b0(1)-b0_wgrid(2);
    b0_wgrid(Nr+1) = 2*b0(Nr)-b0_wgrid(Nr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d2bdz2 = (dbdz(2:Nr+1)-dbdz(1:Nr))/dz; %%% 2nd-order centered difference

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*U.*b0;
    bq2(o,:) = -1i*kx*p0_ugrid*N^2*cosd(topo);
    bq3(o,:) = dpsidz*N^2*sind(topo);
    bq4(o,:) = 1i*kx*p0_ugrid.*dAdz/omega*N^2*sind(topo)*sin(omega*t0);
    bq5(o,:) = kappa*(d2bdz2-kx^2.*b0);

    if(noBQ2)
        bq2(o,:) = 0;
    end
    if(noBQ3)
        bq3(o,:) = 0;
    end
    if(noBQ4)
        bq4(o,:) = 0;
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);

    % for m = 2:Nr
    %     d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    % end
    d2zetadz2(2:Nr) = (z0(1:Nr-1)-2*z0(2:Nr)+z0(3:Nr+1))/dz^2;

    %%%%%%%%%%%% B.C.-11 %%%%%%%%%%%%
    %%% Linear extrapolation
    d2zetadz2 = interp1(zz_wgrid(2:Nr),d2zetadz2(2:Nr),zz_wgrid,'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zq1(o,:) = -1i*kx*U_wgrid.*z0;
    zq2(o,:) = 1i*kx*b0_wgrid*cosd(topo);
    zq3(o,:) = -dbdz*sind(topo);
    zq4(o,:) = nu*(d2zetadz2-kx^2.*z0);

    if(noZQ2)
        zq2(o,:) = 0;
    end
    if(noZQ3)
        zq3(o,:) = 0;
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);



        k_2b = dbdt(o,:);
        k_2z = dzetadt(o,:);

        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o,:)+0.5*dt*k_2b;
        z_3 = zeta(o,:)+0.5*dt*k_2z;
    %%%%%%%%%%%% B.C.-4 %%%%%%%%%%%%
        % z_3(1) = delta/Hshear*cos(t0); %%% No total stress
        z_3(1) = 0; z_3(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
            

    %%%%%%%%%%%% B.C.-6 %%%%%%%%%%%%
    z0(1) = 0; z0(Nr+1) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Solve the boundary value problem using the bvp4c solver.
    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%% Form Initial Guess
        % options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NMax', 20000);
        % p0_guess = real(p0);
        % solinit.y(1,:) = p0_guess;
        % solinit.y(2,2:end-1) = (p0_guess(3:end)-p0_guess(1:end-2))/dz/2; %%% Centered difference
        % solinit.y(2,1)=(p0_guess(2)-p0_guess(1))/dz;
        % solinit.y(2,end)=(p0_guess(end)-p0_guess(end-1))/dz;
    xmesh = zz_wgrid;
    solinit = bvpinit(xmesh, [0 0]);
    options = bvpset('NMax', 10000);
    if(~hydrostatic)
        %%% non-hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit,options);
    else
        %%% hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun_hydro(z,y,zeta0), @bcfun, solinit,options);
    end
    psi0 = sol1.y(1,:);
    zz0 = sol1.x;
    if(length(psi0)>length(p0))
        p0 = interp1(zz0,psi0,zz_wgrid);
        % warning('length(psi0)>length(p0)')
    else
        p0 = psi0;
    end


    %%%%%%%%%%%% B.C.-7 %%%%%%%%%%%%
    %%% Impermeable (w=0) At the ocean bottom
    p0(1) = 0;
    %%% Impermeable (w=0) At the upper boundary 
    p0(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    U = cos(omega*t0)*Atide;
    U_wgrid = cos(omega*t0)*Atide_wgrid;

    % dUdz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz;

    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz; %%% on w-grid

    %%%%%%%%%%%% B.C.-8 %%%%%%%%%%%%
    % %%% No buoyancy flux, Ocean bottom, when diffusivity is not zero
    % %%% dbdz+dB/dz+dB0/dz = 0 ==>  dbdz = -(dB/dz+dB0/dz)
    % dbdz(1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    % %%% No buoyancy flux, Upper boundary, when diffusivity is not zero
    % dbdz(Nr+1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    dbdz(1) = 0;
    dbdz(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    %%%%%%%%%%%% B.C.-9 %%%%%%%%%%%%
    b0_wgrid(1) = 2*b0(1)-b0_wgrid(2);
    b0_wgrid(Nr+1) = 2*b0(Nr)-b0_wgrid(Nr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d2bdz2 = (dbdz(2:Nr+1)-dbdz(1:Nr))/dz; %%% 2nd-order centered difference

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*U.*b0;
    bq2(o,:) = -1i*kx*p0_ugrid*N^2*cosd(topo);
    bq3(o,:) = dpsidz*N^2*sind(topo);
    bq4(o,:) = 1i*kx*p0_ugrid.*dAdz/omega*N^2*sind(topo)*sin(omega*t0);
    bq5(o,:) = kappa*(d2bdz2-kx^2.*b0);

    if(noBQ2)
        bq2(o,:) = 0;
    end
    if(noBQ3)
        bq3(o,:) = 0;
    end
    if(noBQ4)
        bq4(o,:) = 0;
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);

    % for m = 2:Nr
    %     d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    % end
    d2zetadz2(2:Nr) = (z0(1:Nr-1)-2*z0(2:Nr)+z0(3:Nr+1))/dz^2;

    %%%%%%%%%%%% B.C.-11 %%%%%%%%%%%%
    %%% Linear extrapolation
    d2zetadz2 = interp1(zz_wgrid(2:Nr),d2zetadz2(2:Nr),zz_wgrid,'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zq1(o,:) = -1i*kx*U_wgrid.*z0;
    zq2(o,:) = 1i*kx*b0_wgrid*cosd(topo);
    zq3(o,:) = -dbdz*sind(topo);
    zq4(o,:) = nu*(d2zetadz2-kx^2.*z0);

    if(noZQ2)
        zq2(o,:) = 0;
    end
    if(noZQ3)
        zq3(o,:) = 0;
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);






        k_3b = dbdt(o,:);
        k_3z = dzetadt(o,:);


        % Mid-point predictor advancing dt:
        b_4 = buoy(o,:)+dt*k_3b;
        z_4 = zeta(o,:)+dt*k_3z;
    %%%%%%%%%%%% B.C.-5 %%%%%%%%%%%%
        % z_4(1) = delta/Hshear*cos(t0); %%% No total stress
        z_4(1) = 0; z_4(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
         

    %%%%%%%%%%%% B.C.-6 %%%%%%%%%%%%
    z0(1) = 0; z0(Nr+1) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Solve the boundary value problem using the bvp4c solver.
    zeta0 = @(z) interp1(linspace(min(zspan),max(zspan),numel(z0)), z0, z);    
    %%% Form Initial Guess
        % options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-3, 'NMax', 20000);
        % p0_guess = real(p0);
        % solinit.y(1,:) = p0_guess;
        % solinit.y(2,2:end-1) = (p0_guess(3:end)-p0_guess(1:end-2))/dz/2; %%% Centered difference
        % solinit.y(2,1)=(p0_guess(2)-p0_guess(1))/dz;
        % solinit.y(2,end)=(p0_guess(end)-p0_guess(end-1))/dz;
    xmesh = zz_wgrid;
    solinit = bvpinit(xmesh, [0 0]);
    options = bvpset('NMax', 10000);
    if(~hydrostatic)
        %%% non-hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun(z,y,kx,zeta0), @bcfun, solinit,options);
    else
        %%% hydrostctic
        sol1 = bvp4c(@(z,y)bvpfun_hydro(z,y,zeta0), @bcfun, solinit,options);
    end
    psi0 = sol1.y(1,:);
    zz0 = sol1.x;
    if(length(psi0)>length(p0))
        p0 = interp1(zz0,psi0,zz_wgrid);
        % warning('length(psi0)>length(p0)')
    else
        p0 = psi0;
    end


    %%%%%%%%%%%% B.C.-7 %%%%%%%%%%%%
    %%% Impermeable (w=0) At the ocean bottom
    p0(1) = 0;
    %%% Impermeable (w=0) At the upper boundary 
    p0(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    U = cos(omega*t0)*Atide;
    U_wgrid = cos(omega*t0)*Atide_wgrid;

    % dUdz = (U_wgrid(2:Nr+1)-U_wgrid(1:Nr))/dz;

    dbdz(2:Nr) = (b0(2:Nr)-b0(1:Nr-1))/dz; %%% on w-grid

    %%%%%%%%%%%% B.C.-8 %%%%%%%%%%%%
    % %%% No buoyancy flux, Ocean bottom, when diffusivity is not zero
    % %%% dbdz+dB/dz+dB0/dz = 0 ==>  dbdz = -(dB/dz+dB0/dz)
    % dbdz(1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    % %%% No buoyancy flux, Upper boundary, when diffusivity is not zero
    % dbdz(Nr+1) = - (omega/Shear) * (delta/Hshear) * (cosd(topo)/sind(topo)) ...
    %     +(delta/Hshear)*sin(t0); 
    dbdz(1) = 0;
    dbdz(Nr+1) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b0_wgrid(2:Nr) = 0.5*(b0(2:Nr)+b0(1:Nr-1));

    %%%%%%%%%%%% B.C.-9 %%%%%%%%%%%%
    b0_wgrid(1) = 2*b0(1)-b0_wgrid(2);
    b0_wgrid(Nr+1) = 2*b0(Nr)-b0_wgrid(Nr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d2bdz2 = (dbdz(2:Nr+1)-dbdz(1:Nr))/dz; %%% 2nd-order centered difference

    p0_ugrid = 0.5*(p0(1:Nr)+p0(2:Nr+1));
    dpsidz = (p0(2:Nr+1)-p0(1:Nr))/dz;

    bq1(o,:) = -1i*kx*U.*b0;
    bq2(o,:) = -1i*kx*p0_ugrid*N^2*cosd(topo);
    bq3(o,:) = dpsidz*N^2*sind(topo);
    bq4(o,:) = 1i*kx*p0_ugrid.*dAdz/omega*N^2*sind(topo)*sin(omega*t0);
    bq5(o,:) = kappa*(d2bdz2-kx^2.*b0);

    if(noBQ2)
        bq2(o,:) = 0;
    end
    if(noBQ3)
        bq3(o,:) = 0;
    end
    if(noBQ4)
        bq4(o,:) = 0;
    end

    dbdt(o,:) = bq1(o,:) + bq2(o,:) + bq3(o,:) ...
              + bq4(o,:) + bq5(o,:);

    % for m = 2:Nr
    %     d2zetadz2(m) = (z0(m-1)-2*z0(m)+z0(m+1))/dz^2;
    % end
    d2zetadz2(2:Nr) = (z0(1:Nr-1)-2*z0(2:Nr)+z0(3:Nr+1))/dz^2;

    %%%%%%%%%%%% B.C.-11 %%%%%%%%%%%%
    %%% Linear extrapolation
    d2zetadz2 = interp1(zz_wgrid(2:Nr),d2zetadz2(2:Nr),zz_wgrid,'linear','extrap');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zq1(o,:) = -1i*kx*U_wgrid.*z0;
    zq2(o,:) = 1i*kx*b0_wgrid*cosd(topo);
    zq3(o,:) = -dbdz*sind(topo);
    zq4(o,:) = nu*(d2zetadz2-kx^2.*z0);

    if(noZQ2)
        zq2(o,:) = 0;
    end
    if(noZQ3)
        zq3(o,:) = 0;
    end

    dzetadt(o,:) = zq1(o,:) + zq2(o,:) + zq3(o,:) + zq4(o,:);



        k_4b = dbdt(o,:);
        k_4z = dzetadt(o,:);
        % Simpson rule corrector advancing dt:
        buoy(o+1,:) = buoy(o,:) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1,:) = zeta(o,:) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;


    elseif(useAB3)
        %%% Third-order Adams-Bashforth method %%%
        if (o <= 2)
            %%% Use RK4 for the first 2 time steps
            RK4;
            %%% Use Euler-forward for the first 2 time steps
            % buoy(o+1,:) = buoy(o,:) + dbdt(o,:)*dt;
            % zeta(o+1,:) = zeta(o,:) + dzetadt(o,:)*dt;
        else
            buoy(o+1,:) = buoy(o,:) + dt*( (23/12)*dbdt(o,:)    - (16/12)*dbdt(o-1,:)    + (5/12)*dbdt(o-2,:) );
            zeta(o+1,:) = zeta(o,:) + dt*( (23/12)*dzetadt(o,:) - (16/12)*dzetadt(o-1,:) + (5/12)*dzetadt(o-2,:) );
        end
    else
        %%% Euler forward %%%
        buoy(o+1,:) = dbdt(o,:)*dt;
        zeta(o+1,:) = dzetadt(o,:)*dt;
    end

    %%%%%%%%%%%% B.C.-2 %%%%%%%%%%%%
    %%% No-stress (free-slip)
    zeta(o+1,1) = 0; zeta(o+1,Nr+1) = 0; 

    %%% No total stress at the ocean bottom
    % zeta(o,1) = cos(t0*omega);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(rem(o,round(Nt/10))==0)
        Progress = o/Nt
        % outputname_mid = [expdir 'output_' num2str(round(Progress*10)) '.mat'];
        % re_psi = real(psi);  
        % re_zeta = real(zeta);
        % re_buoy = real(buoy); 
        % uuu = -real((psi(:,2:Nr+1)-psi(:,1:Nr))/dz);
        % www = real(1i*kx*psi);
        % save(outputname_mid,'re_psi','re_zeta','re_buoy','uuu','www','zz','tt')
        % clear re_psi re_zeta re_buoy uuu www
    end

end

re_psi = real(psi);   re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;
re_dbdz = real(dbdz);

re_bq1 = real(bq1);
re_bq2 = real(bq2);
re_bq3 = real(bq3);
re_bq4 = real(bq4);
re_bq5 = real(bq5);

re_zq1 = real(zq1);
re_zq2 = real(zq2);
re_zq3 = real(zq3);
re_zq4 = real(zq4);


dbuoydz = zeros(Nt,Nr);
for m = 2:Nr-1
    dbuoydz(:,m) = (re_buoy(:,m+1)-re_buoy(:,m-1))/dz;
end


plot_tidx = 1:10:Nt;

% load_colors;

DIV = 1;
uuu = -real((psi(:,2:Nr+1)-psi(:,1:Nr))/dz);
www = real(1i*kx*psi);

% h=figure(5);
% set(h,'color','w','Visible', FigureIsVisible,'Position',[67 346 1015 619]);
% clf;
% subplot(3,2,1)
% pcolor(tt(plot_tidx)/t1hour,zz_wgrid,re_psi(plot_tidx,:)');shading flat;colorbar;
% colormap(redblue)
% set(gca,'Fontsize',fontsize);
% ylabel('HAB (m)');
% title('Streamfunction \psi','Fontsize',fontsize+3);
% set(gca,'color',gray);
% aaa = max(max(abs(re_psi))/DIV);
% if(~isnan(aaa))
% caxis([-1 1]*aaa)
% end
% 
% subplot(3,2,2)
% pcolor(tt(plot_tidx)/t1hour,zz_wgrid,re_zeta(plot_tidx,:)');shading flat;colorbar;
% set(gca,'Fontsize',fontsize);
% ylabel('HAB (m)');
% title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
% set(gca,'color',gray);
% aaa = max(max(abs(re_zeta))/DIV);
% if(~isnan(aaa))
% caxis([-1 1]*aaa)
% end
% 
% subplot(3,2,3)
% pcolor(tt(plot_tidx)/t1hour,zz,re_buoy(plot_tidx,:)');shading flat;colorbar;
% set(gca,'Fontsize',fontsize);
% ylabel('HAB (m)');
% title('Buoyancy perturbation b^\prime','Fontsize',fontsize+3);
% set(gca,'color',gray);
% aaa = max(max(abs(re_buoy))/DIV);
% if(~isnan(aaa))
% caxis([-1 1]*aaa)
% end
% 
% 
% subplot(3,2,5)
% pcolor(tt(plot_tidx)/t1hour,zz,uuu(plot_tidx,:)');shading flat;colorbar;
% colormap(redblue)
% set(gca,'Fontsize',fontsize);
% ylabel('HAB (m)');xlabel('Time (hours)')
% title('Horizontal velocity perturbation u^\prime','Fontsize',fontsize+3);
% set(gca,'color',gray);
% aaa = max(max((abs(uuu))))/DIV;
% if(~isnan(aaa)&& aaa~=0)
% caxis([-1 1]*aaa)
% end
% 
% subplot(3,2,6)
% pcolor(tt(plot_tidx)/t1hour,zz_wgrid,www(plot_tidx,:)');shading flat;colorbar;
% set(gca,'Fontsize',fontsize);
% ylabel('HAB (m)');xlabel('Time (hours)')
% title('Vertical velocity w','Fontsize',fontsize+3);
% set(gca,'color',gray);
% if(kx~=0)
%     caxis([-1 1]*max(max((abs(www))))/DIV)
% end
% s = struct('h',h);
% save(sprintf([expdir 'fig5.png']),"-fromstruct",h);

zidx = 1:Nr;
% zidx = 63:Nr-62;
% zidx = 125:Nr;
% zidx = 250:Nr;

bq1_int = real((sum(bq1(:,zidx),2)))';
bq2_int = real((sum(bq2(:,zidx),2)))';
bq3_int = real((sum(bq3(:,zidx),2)))';
bq4_int = real((sum(bq4(:,zidx),2)))';
bq5_int = real((sum(bq5(:,zidx),2)))';

% bq_all = bq1_int+bq2_int+bq3_int+bq4_int+bq5_int;
bq_all = 1;

zidx = 1:Nr+1;

zq1_int = real((sum(zq1(:,zidx),2)))';
zq2_int = real((sum(zq2(:,zidx),2)))';
zq3_int = real((sum(zq3(:,zidx),2)))';
zq4_int = real((sum(zq4(:,zidx),2)))';
% zq_all = zq1_int+zq2_int+zq3_int+zq4_int;
zq_all = 1;

lw = 2;
% close all

fit_span = round(Nt/NTtide*3)+1:Nt-1;
% fit_span = 1:Nt;

% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TPE = 0;
KE_PE = TKE+TPE;

KE_PE_zavg = mean(KE_PE,2)';
xxplot = tt/t1hour;
yyplot = log(KE_PE_zavg)/2;
[pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% GrowthRate_KE(Nexp_lambda,Nexp_shear) = pKE(1);
pKE(1);
[y_fit,delta_fit] = polyval(pKE,xxplot,S);

b2 = mean(re_buoy.^2,2)';
yyplot_b2 = log(b2)/2;
[pb2,S_b2] = polyfit(xxplot(fit_span),yyplot_b2(fit_span),1); 
% GrowthRate(Nexp_lambda,Nexp_shear) = pb2(1);
[y_fit_b2,delta_fit_b2] = polyval(pb2,xxplot,S_b2);


% h=figure(8);
% clf;
% set(h,'color','w','Visible', FigureIsVisible,'Position',[85 222 979 420]);
% plot(xxplot/12,yyplot,'LineWidth',2)
% hold on
% plot(xxplot/12,yyplot_b2,'LineWidth',2)
% plot(xxplot(fit_span)/12,y_fit(fit_span),':','LineWidth',2)
% plot(xxplot(fit_span)/12,y_fit_b2(fit_span),':','LineWidth',2)
% grid on;grid minor;
% set(gca,'Fontsize',20);
% % ylim([pKE(2)-3 pKE(2)+pKE(1)*max(xxplot)+2])
% xlabel('$t$ (tidal cycle)','Interpreter','Latex')
% ylabel('$\ln(e)/2$','Interpreter','Latex')
% hold off;axis tight
% legend('TKE','(b^\prime)^2','Position',[0.8141 0.1988 0.0684 0.1393])
% s = struct('h',h);
% save(sprintf([expdir 'KE.png']),"-fromstruct",h);

% close all;
% save(outputname)

    % clear b0 b_wgrid b0_wgrid b_2 b_3 b_4 p0 p0_ugrid psi psi0 sol1 solinit ...
    %     z0 z_2 z_3 z_4 zeta dbdz ...
    %     d2bdz2 d2psidz2 d2zetadz2  dpsidz dUtidedz dzetadt ...
    %     bq1 bq2 bq3 bq4 bq5 zq1 zq2 zq3 zq4 ...
    %     k_1b k_1z k_2b k_2z k_3b k_3z k_4b k_4z h ...
    %     buoy re_zq1 re_zq2 re_zq3 re_zq4 re_bq1 re_bq2 re_bq3 re_bq4 re_bq5 Utide ...
    %     uuu  re_dbdz re_d2bdz2 re_d2zetadz2 
    
    % save(outputname)
    s = struct('buoy',buoy,'zeta',zeta,'psi',psi, ...
            're_buoy',re_buoy,'re_psi',re_psi,'pKE',pKE,'pb2',pb2,...
            'uuu',uuu,'www',www,'NTtide',NTtide,'Nr',Nr,'Nt',Nt,'Utide',Utide,...
            'tt',tt,'t1hour',t1hour,'zz',zz,'dz',dz,'nu',nu,'kappa',kappa,'fit_span',fit_span,...
            'bq1_int',bq1_int,'bq2_int',bq2_int,'bq3_int',bq3_int,'bq4_int',bq4_int,'bq5_int',bq5_int,...
            'zq1_int',zq1_int,'zq2_int',zq2_int,'zq3_int',zq3_int,'zq4_int',zq4_int,...
            'Shear',Shear,'lambda',lambda,'topo',topo,'Atide',Atide,'dAdz',dAdz...
            );
            
    save(sprintf(outputname),"-fromstruct",s);


    end
    
end


    %%% Code equations: see https://www.mathworks.com/help/matlab/math/solve-bvp-with-two-solutions.html
    function dydz = bvpfun(z,y,kx,zeta)
        dydz = [y(2)
            kx^2*y(1)+zeta(z)];
    end

    function dydz = bvpfun_hydro(z,y,zeta)
        dydz = [y(2)
              zeta(z)];
    end

    %%%%%%%%%%%% B.C.-12 %%%%%%%%%%%%
    %%% Code boundary conditions for the streamfunction: phi = 0 at z=0 and z=1
    function res = bcfun(ya,yb)
        res = [ya(1)
               yb(1)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

