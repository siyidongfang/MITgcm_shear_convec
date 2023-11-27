
clear;
% close all;

addpath ../analysis/colormaps/

useRK4 = true;

%%%% Define constants
Hdepth = 1500;
Hshear = 300;
Shear = 1*0.8e-3; 
N = 1e-3;
topo = 4;
omega = 2*pi/43200;
nu = 2e-4;
kappa = 2e-4;
Pr = nu/kappa;
m1km = 1000;
t1hour = 3600;

delta = Hdepth;
% delta = sqrt(2*nu/omega);
% delta = U0/omega;

%%% Model dimension
Lz = Hdepth/delta;     % dimensionless domain height
% dz = 10;               % dimensionless vertical grid spacing
dz = 0.01;
Nr = round(Lz/dz)+2;
% zz = dz/2:dz:(Nr*dz-dz/2);  % Height above topography
zz = 0:dz:((Nr-1)*dz);

Lshear = Hshear/delta; % dimensionless vertical scale for velocity shear
Nshear = round(Lshear/dz);
Hshear = zz(Nshear)*delta;
U0 = Hshear * Shear;
Re = U0*delta/nu;

NTtide = 3;
Lt = NTtide*43200*omega; % dimensionless simulation time
dt = 0.01;
Nt = round(Lt/dt);
tt = dt:dt:Nt*dt;

C1 = U0/delta/omega;
C2 = N*sind(topo)/omega;
C3 = nu/delta^2/omega;
C4 = kappa/delta^2/omega;


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
d2zetadz2 = zeros(1,Nr);
dUdz = zeros(1,Nr);
U = zeros(1,Nr);

lambda = 1*m1km;
kx = 2*pi/(lambda/delta) %%% Dimensionless wave number kx = 0.0000001, no unstable layer

%%% Initial condition
% buoy(1,:) = 0;
% buoy(1,:) = (rand(1,Nr)-1/2)/1.79e300;
buoy(1,:) = (rand(1,Nr)-1/2)/1e20;
% buoy(1,:) = 1/1e20;
psi(1,:) = 0;
zeta(1,:) = 0;
% psi(1,:) = (rand(1,Nr)-1/2)/1.79e30;
% for m = 2:Nr-1
%     dpsidz(m)   = (psi(1,m+1)-psi(1,m-1))/2/dz;
%     d2psidz2(m) = (psi(1,m-1)-2*psi(1,m)+psi(1,m+1))/dz^2;
% end
% zeta(1,:) = d2psidz2-kx^2*psi(1,:);

%%% Background tidal velocity (dimensionless)
Atide = zeros(1,Nr);
for m=1:Nshear
    Atide(m) = Shear*zz(m)*delta;
end
for m=Nshear+1:Nr
    Atide(m) = U0;
end
% idx_smooth = Nshear-round(Nshear/4):Nshear+4*round(Nshear/4);
idx_smooth = 1:Nr;
Atide(idx_smooth) = smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(Atide(idx_smooth)))))))))))))))))))))';
% Atide(idx_smooth) = smooth(smooth(smooth(Atide(idx_smooth))))';

Utide =cos(tt)'.*Atide/U0;

dUtidedz = zeros(Nt,Nr);
for m = 2:Nr-1
    dUtidedz(:,m)   = (Utide(:,m+1)-Utide(:,m-1))/2/dz;
end



% figure(1);clf;
% plot(Atide,zz*delta)
% 
% figure(2);clf;
% plot(diff(Atide)./diff(zz)/delta,0.5*(zz(1:end-1)+zz(2:end))*delta)
% grid on;grid minor;
% 
% figure(3);clf;
% pcolor(tt/omega/3600,zz*delta,Utide'*U0);shading flat;colormap redblue; colorbar;
% 
% figure(4);clf;
% pcolor(tt/omega/3600,zz*delta,dUtidedz'*U0/delta);shading flat;colormap redblue; colorbar;

%%

%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Fourth-order Runge-Kutta or Euler forward scheme for time advancement

% buoy(1,:) = rand(1,Nr)/1e20;
% buoy(1,:) = -[1:Nr]/1e20;

p0 = psi(1,:);

for o=1:Nt-1

    t0 = tt(o);
    b0 = buoy(o,:);
    z0 = zeta(o,:);
    
    tendency;
    psi(o,:) = p0;

    if(useRK4)
    %%% Fourth-order Runge-Kutta %%%
        k_1b = dbdt(o,:);
        k_1z = dzetadt(o,:);
        % Euler forward predictor advancing dt/2:
        b_2 = buoy(o,:)+0.5*dt*k_1b;
        z_2 = zeta(o,:)+0.5*dt*k_1z;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        z0 = z_2;
        tendency;
        k_2b = dbdt(o,:);
        k_2z = dzetadt(o,:);
        % Euler backward corrector advancing dt/2:
        b_3 = buoy(o,:)+0.5*dt*k_2b;
        z_3 = zeta(o,:)+0.5*dt*k_2z;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        z0 = z_3;
        tendency;
        k_3b = dbdt(o,:);
        k_3z = dzetadt(o,:);
        % Mid-point predictor advancing dt:
        b_4 = buoy(o,:)+dt*k_3b;
        z_4 = zeta(o,:)+dt*k_3z;
        t0 = tt(o)+dt;
        b0 = b_4;
        z0 = z_4;
        tendency;
        k_4b = dbdt(o,:);
        k_4z = dzetadt(o,:);
        % Simpson rule corrector advancing dt:
        buoy(o+1,:) = buoy(o,:) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        zeta(o+1,:) = zeta(o,:) + (1/6)*(k_1z+2*k_2z+2*k_3z+k_4z)*dt;
    else
    %%% Euler forward %%%
        buoy(o+1,:) = dbdt(o,:)*dt;
        zeta(o+1,:) = dzetadt(o,:)*dt;
    end
end

re_psi = real(psi); re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;


%%%%% Floquet stability
oT = round(43200*omega/dt);% The time step after one tidal cycle

muk_psi = re_psi(oT*2-1,:)./re_psi(oT-1,:);
muk_zeta = re_zeta(oT*2-1,:)./re_zeta(oT-1,:);
muk_buoy = re_buoy(oT*2-1,:)./re_buoy(oT-1,:);

sum(muk_psi>1)
sum(muk_zeta>1)
sum(muk_buoy>1)

%%% Dimensional veriables
zzd = zz*delta;    
ttd = tt/omega;     
re_psid = re_psi*U0/delta;
re_zetad = re_zeta*U0*delta;
re_buoyd = re_buoy*N^2*sind(topo)/omega;

dbuoydz = zeros(Nt,Nr);
for m = 2:Nr-1
    dbuoydz(:,m) = (re_buoyd(:,m+1)-re_buoyd(:,m-1))/2/dz/delta;
end

%%
plot_tidx = 1:1:Nt;
load_colors;


figure(5)
fontsize = 16;
clf;set(gcf,'color','w','Position',[44 241 654 728]);
subplot(3,1,1)
pcolor(ttd(plot_tidx)/t1hour,zzd,re_psid(plot_tidx,:)');shading flat;colorbar;
colormap(redblue)
% colormap(WhiteBlueGreenYellowRed(0));
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Streamfunction \psi','Fontsize',fontsize+3);
set(gca,'color',gray);
% clim([-1 1]/1e10)
clim([-1 1]*U0*delta*delta)
% aaa = max(max(abs(re_psid)));
% clim([-1 1]*aaa/100)

subplot(3,1,2)
pcolor(ttd(plot_tidx)/t1hour,zzd,re_zetad(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
set(gca,'color',gray);
% clim([-1 1]/1e6)
clim([-1 1]*U0*delta)
% aaa = max(max(abs(re_zetad)));
% clim([-1 1]*aaa/100)

subplot(3,1,3)
pcolor(ttd(plot_tidx)/t1hour,zzd,dbuoydz(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Stratification perturbation','Fontsize',fontsize+3);
set(gca,'color',gray);
% clim([-1 1]/1e10)
clim([-1 1]*N^2*sind(topo)/omega)
% aaa = max(max(abs(re_buoyd)));
% clim([-1 1]*aaa/100)

fname = ['exps_instability/Shear' num2str(Shear) '_lambda' num2str(lambda) '.mat'];
save(fname,'buoy','zeta','psi', ...
    'dbdt','dzetadt', ...
    'bq1','bq2','bq3','bq4','bq5','bq6',...
    'zq1','zq2','zq3','zq4','zq5')




