
clear;close all;

addpath ../analysis/colormaps/

%%%% Define constants
Hdepth = 1500;
Hshear = 300;
Shear = 1*0.8e-3; 
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
dz = 5;               % dimensionless vertical grid spacing
Nr = round(Lz/dz)+1;
zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography

Lshear = Hshear/delta; % dimensionless vertical scale for velocity shear
Nshear = round(Lshear/dz);
Hshear = zz(Nshear)*delta;
U0 = Hshear * Shear;
Re = U0*delta/nu;
D = Re*delta/2/Hshear;

Lt = 0.5*86400*omega; % dimensionless simulation time
dt = 0.001;
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


%%
%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Euler forward scheme for time advancement
% %%%%  TO DO: Fourth-order Runge-Kutta scheme for time advancement

%%%% Ignore diffusion

for o=1:Nt-1

    psi(o,:) = cumsum(cumsum(zeta(o,:)*dz)*dz)-Utide(o,1)*Lz;

    for m = 2:Nr-1
        dpsidz(m)   = (psi(o,m+1)-psi(o,m-1))/2/dz;
        d2psidz2(m) = (psi(o,m-1)-2*psi(o,m)+psi(o,m+1))/dz^2;
        dbdz(m)   = (buoy(o,m+1)-buoy(o,m-1))/2/dz;
        d2bdz2(m) = (buoy(o,m-1)-2*buoy(o,m)+buoy(o,m+1))/dz^2;
    end

    dzetadt = 1i*kx*Re/2*dUtidedz(o,:).*dpsidz ...
         -1i*kx*Re/2*Utide(o,:).*(d2psidz2-kx^2*psi(o,:)) ...
         +C^2*(1i*kx*cotd(topo)-dbdz);

    dbdt = dpsidz -1i*kx*cotd(topo)*psi(o,:)...
         +1i*kx*Re/2*dUtidedz(o,:).*tan(o*dt).*psi(o,:) ...
         -1i*kx*Re/2*Utide(o,:).*buoy(o,:);

    % dzetadt(zidx) = 1i*kx*D*(dpsidz(zidx)-zz(zidx).*(d2psidz2(zidx)-kx^2*psi(o,zidx)))*cos(o*dt) ...
    %      +C^2*(1i*kx*cotd(topo)-dbdz(zidx));
    % 
    % dbdt(zidx) = dpsidz(zidx) + (-1i*kx*cotd(topo)+1i*kx*D*sin(o*dt))*psi(o,zidx) ...
    %      -1i*kx*D*zz(zidx)*cos(o*dt).*buoy(o,zidx);

    %%% Consider dissipation
    dbdt = dbdt + (d2bdz2 - kx^2.*buoy(o,:))/2/Pr;

    %%% Boundary condition
    psi(o,end)=0;

    %%% Euler forward
    buoy(o+1,:) = dbdt*dt;
    zeta(o+1,:) = dzetadt*dt;

    %%% Fourth-order Runge-Kutta


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

subplot(3,1,2)
pcolor(ttd/3600,zzd,re_zetad');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
clim([-4 4]/1e4*U0*delta)


subplot(3,1,3)
pcolor(ttd/3600,zzd,re_buoyd');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Buoyancy perturbation b','Fontsize',fontsize+3);
clim([-1000 1000]*N^2*sind(topo)/omega)








