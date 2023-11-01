
%%% Compare three terms in the temperature equation:
%%% u*dT/dx, v*dT/dy(==0), w*dT/dz, w*N^2*cos(theta)/g/alpha, u*N^2*sin(theta)/g/alpha

%%%% TO DO: use UVELTH, WVELTH??
%%%% TO DO: use phase-average??
%%%% TO DO: plot average over a tidal cycle??
%%%% TO DO: use hourly-mean??
%%%% TO DO: plot temperature tendency, diffusion, and the total??


clear;close all;
ne =1;
load_all

topo_slope = 4;
cos_slope = cosd(topo_slope);
sin_slope = sind(topo_slope);
N2_const = 1e-6;

CLIM = [-1 1]/5e4;
YLIM = [1100 1500];XLIM = [-Lx/2/1000 Lx/2/1000];


Hz = sum(delR);
N2const = (1e-3)^2;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);

for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end

o1 = 36;
o2 = 240;

for o=o1:o2
nIter = dumpIters(o);
time_h = nIter.*deltaT./3600;

tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
uu = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
vv = squeeze(rdmds([exppath,'/results/VVEL_inst'],nIter));
ww = squeeze(rdmds([exppath,'/results/WVEL_inst'],nIter));

wN = ww*N2_const*cos_slope/gravity/tAlpha;
uN = uu*N2_const*sin_slope/gravity/tAlpha;

dx = delX(1);
dtdx = (tt([Nx 1:Nx-1],:) - tt)/dx; %%% on u-grid
udtdx = uu.*dtdx;

dtdz = zeros(Nx,Nr);
dtdz(:,2:Nr) = -diff(tt,1,2)./(-diff(zz)); %%% on w-grid
wdtdz = ww.*dtdz;


tt = tt + tt_background;
rho = rhoConst.*(1-(tt-tRef)*tAlpha);
N2 = NaN*zeros(Nx,Nr);
N2(:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));


%%% Plot the 4 terms
figure(1)
clf;set(gcf,'color','w','Position',[101 63 1051 788]);
subplot(3,2,1)
pcolor(xx/1000,-zz,udtdx');colormap(redblue);
hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
title(['$u\,\partial_x\theta$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
clim(CLIM);ylim(YLIM);xlim(XLIM);

subplot(3,2,2)
pcolor(xx/1000,-zz,wdtdz');colormap(redblue);
hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
title(['$w\,\partial_z\theta$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
clim(CLIM);ylim(YLIM);xlim(XLIM);

subplot(3,2,3)
pcolor(xx/1000,-zz,uN');colormap(redblue);
hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
title(['$u N^2\sin\theta/(g\alpha)$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
clim(CLIM);ylim(YLIM);xlim(XLIM);

subplot(3,2,4)
pcolor(xx/1000,-zz,wN');colormap(redblue);
hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
title(['$w N^2\cos\theta/(g\alpha)$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
clim(CLIM);ylim(YLIM);xlim(XLIM);

% subplot(3,2,5) %%% Temperature tendency
% subplot(3,2,6) %%% Diffusion

subplot(3,2,5)
pcolor(xx/1000,-zz,N2')
shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
title('$N^2\ (s^{-2})$','Fontsize',fontsize+3,'interpreter','latex')
clim([-2 2]/1e6)
ylim(YLIM)

end




