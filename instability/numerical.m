
clear;close all;

addpath ../analysis/colormaps/

useRK4 = false;

%%%% Define constants
Hdepth = 1500;
Hshear = 300;
Shear = 0.4*0.8e-3; 
N = 1e-3;
topo = 4;
omega = 2*pi/43200;
nu = 2e-6;
kappa = 2e-6;
Pr = nu/kappa;
delta = sqrt(2*nu/omega);
C = N*sind(topo)/omega;

%%% Model dimension
Lz = Hdepth/delta;     % dimensionless domain height
dz = 100;               % dimensionless vertical grid spacing
Nr = round(Lz/dz)+1;
zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography

Lshear = Hshear/delta; % dimensionless vertical scale for velocity shear
Nshear = round(Lshear/dz);
Hshear = zz(Nshear)*delta;
U0 = Hshear * Shear;
Re = U0*delta/nu;

Lt = 2*43200*omega; % dimensionless simulation time
dt = 0.0005;
Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;

%%%% Define variables
psi = zeros(Nt,Nr);
zeta = zeros(Nt,Nr);
buoy = zeros(Nt,Nr);

dbdt = zeros(1,Nr);
dzetadt = zeros(1,Nr);

dbdz = zeros(1,Nr);
d2bdz2 = zeros(1,Nr);
dpsidz = zeros(1,Nr);
d2psidz2 = zeros(1,Nr);
dUdz = zeros(1,Nr);
U = zeros(1,Nr);

kx = 0.01; %%% wave number

%%% Initial condition
psi(1,:) = 0;
zeta(1,:) = 0;
buoy(1,:) = 0;


%%% Background tidal velocity (dimensionless)
Atide = zeros(1,Nr);
for m=1:Nshear
    Atide(m) = Shear*zz(m)*delta;
end
for m=Nshear+1:Nr
    Atide(m) = U0;
end
idx_smooth = Nshear-round(Nshear/4):Nshear+4*round(Nshear/4);
Atide(idx_smooth) = smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(Atide(idx_smooth)))))))))))))))))))))';

Utide =cos(tt)'.*Atide/U0;

dUtidedz = zeros(Nt,Nr);
for m = 2:Nr-1
    dUtidedz(:,m)   = (Utide(:,m+1)-Utide(:,m-1))/2/dz;
end

figure(1);clf;
plot(Atide,zz*delta)

figure(2);clf;
plot(diff(Atide)./diff(zz)/delta,0.5*(zz(1:end-1)+zz(2:end))*delta)

figure(3);clf;
pcolor(tt/omega/3600,zz*delta,Utide'*U0);shading flat;colormap redblue; colorbar;

figure(4);clf;
pcolor(tt/omega/3600,zz*delta,dUtidedz'*U0/delta);shading flat;colormap redblue; colorbar;



%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Euler forward scheme for time advancement
%%%% Fourth-order Runge-Kutta scheme for time advancement

%%%% Ignore diffusion

for o=1:Nt-1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fourth-order Runge-Kutta %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t0 = tt(o);
    U = cos(t0)*Atide/U0;
    psi(o,:) = cumsum(cumsum(zeta(o,:)*dz)*dz)-U(1)*Lz;
    psi(o,end)=0;     % Boundary conditions
    psi(o,1)=0;
    buoy(o,1) = buoy(o,2);
    buoy(o,Nr) = buoy(o,Nr-1);

    p0 = psi(o,:);
    b0 = buoy(o,:);
    z0 = zeta(o,:);

    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);
    k_1b = dbdt;
    k_1z = dzetadt;

    if(useRK4)
    %%% Euler forward predictor advancing dt/2:
    b_2 = buoy(o,:)+0.5*dt*k_1b;
    z_2 = zeta(o,:)+0.5*dt*k_1z;

    t0 = tt(o)+dt/2;
    b0 = b_2;
    z0 = z_2;

    U = cos(t0)*Atide/U0;
    p0 = cumsum(cumsum(z0*dz)*dz)-U(1)*Lz;
    p0(1)=0;
    p0(Nr)=0;     
    b0(1) = b0(2); 
    b0(Nr) = b0(Nr-1); 
    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);

    k_2b = dbdt;
    k_2z = dzetadt;

    %%% Euler backward corrector advancing dt/2:
    b_3 = buoy(o,:)+0.5*dt*k_2b;
    z_3 = zeta(o,:)+0.5*dt*k_2z;

    t0 = tt(o)+dt/2;
    b0 = b_3;
    z0 = z_3;
    U = cos(t0)*Atide/U0;
    p0 = cumsum(cumsum(z0*dz)*dz)-U(1)*Lz;
    p0(1)=0;
    p0(Nr)=0;     
    b0(1) = b0(2); 
    b0(Nr) = b0(Nr-1); 
    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);

    k_3b = dbdt;
    k_3z = dzetadt;


    %%% Mid-point predictor advancing dt:
    b_4 = buoy(o,:)+dt*k_3b;
    z_4 = zeta(o,:)+dt*k_3z;

    t0 = tt(o)+dt;
    b0 = b_4;
    z0 = z_4;
    U = cos(t0)*Atide/U0;
    p0 = cumsum(cumsum(z0*dz)*dz)-U(1)*Lz;
    p0(1)=0;
    p0(Nr)=0;     
    b0(1) = b0(2); 
    b0(Nr) = b0(Nr-1); 
    for m = 2:Nr-1
        dUdz(m) = (U(m+1)-U(m-1))/2/dz;
        dpsidz(m)   = (p0(m+1)-p0(m-1))/2/dz;
        d2psidz2(m) = (p0(m-1)-2*p0(m)+p0(m+1))/dz^2;
        dbdz(m)   = (b0(m+1)-b0(m-1))/2/dz;
        d2bdz2(m) = (b0(m-1)-2*b0(m)+b0(m+1))/dz^2;
    end
    dbdt = dpsidz -1i*kx*cotd(topo)*p0...
         +1i*kx*Re/2*dUdz.*tan(t0).*p0 ...
         -1i*kx*Re/2*U.*b0 ...
         +(d2bdz2 - kx^2.*b0)/2/Pr; %%% dissipation
    dzetadt = 1i*kx*Re/2*dUdz.*dpsidz ...
         -1i*kx*Re/2*U.*(d2psidz2-kx^2*p0) ...
         +C^2*(1i*kx*cotd(topo)*b0-dbdz);

    k_4b = dbdt;
    k_4z = dzetadt;

    %%% Simpson rule corrector advancing dt:
    buoy(o+1,:) = buoy(o,:) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
    zeta(o+1,:) = zeta(o,:) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;
    
    else
        %%% Euler forward
        buoy(o+1,:) = dbdt*dt;
        zeta(o+1,:) = dzetadt*dt;
    end

end

re_psi = real(psi); re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;


% %%%%% Floquet stability
% oT = round(43200*omega/dt);% The time step after one tidal cycle
% 
% muk_psi = re_psi(oT*2-100,:)./re_psi(oT-100,:);
% muk_zeta = re_zeta(oT*2-100,:)./re_zeta(oT-100,:);
% muk_buoy = re_buoy(oT*2-100,:)./re_buoy(oT-100,:);
% 
% sum(muk_psi>1)
% sum(muk_zeta>1)
% sum(muk_buoy>1)

%%% Dimensional veriables
zzd = zz*delta;    
ttd = tt/omega;     
re_psid = re_psi*U0/delta;
re_zetad = re_zeta*U0*delta;
re_buoyd = re_buoy*N^2*sind(topo)/omega;

%%

figure(5)
fontsize = 16;
clf;set(gcf,'color','w');
subplot(3,1,1)
pcolor(ttd/3600,zzd,re_psid');shading flat;colorbar;colormap(redblue);
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Streamfunction \psi','Fontsize',fontsize+3);
clim([-100 100]*U0/delta)
% aaa = max(max(abs(re_psid)));
% clim([-1 1]*aaa/1e60)

subplot(3,1,2)
pcolor(ttd/3600,zzd,re_zetad');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
clim([-4 4]/1e4*U0*delta)
% aaa = max(max(abs(re_zetad)));
% clim([-1 1]*aaa/1e60)

subplot(3,1,3)
pcolor(ttd/3600,zzd,re_buoyd');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Buoyancy perturbation b','Fontsize',fontsize+3);
clim([-0.1 0.1]*N^2*sind(topo)/omega)
% aaa = max(max(abs(re_buoyd)));
% clim([-1 1]*aaa/1e60)







