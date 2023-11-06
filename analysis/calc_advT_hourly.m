
%%% Compare three terms in the temperature equation:
%%% u*dT/dx, v*dT/dy(==0), w*dT/dz, w*N^2*cos(theta)/g/alpha, u*N^2*sin(theta)/g/alpha

%%% Use ADVr_TH and ADVx_TH, plot phase average over 10 tidal cycles

%%%% TO DO: plot temperature tendency, diffusion, and the total??


clear;close all;
ne =1;
load_all

%%% Grid spacing matrices
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

DZC = zeros(Nx,Ny,Nr);
DZC(:,:,2:Nr) = repmat(reshape(-diff(zz),[1 1 Nr-1]),[Nx Ny 1]);
Cell_Vol = delX(1).*delY(1).*DZC;
Cell_Vol = squeeze(Cell_Vol);

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

dx = delX(1);

o1 = 1;
o2 = 120;
No = o2-o1+1

%%% Integrate over the bottom 250m
Nbot = round(250/delR(end));
dutdx_int = size(1,No);
dwtdz_int = size(1,No);
tt_adv_int = size(1,No);
tt_tend_int = size(1,No);
uN_int = size(1,No);
wN_int = size(1,No);

for o=o1:o2
    h_tides = o;
    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    
    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    ADVx_TH = squeeze(rdmds([exppath,'/results/ADVx_TH'],nIter));
    ADVr_TH = squeeze(rdmds([exppath,'/results/ADVr_TH'],nIter));
    TOTTTEND = squeeze(rdmds([exppath,'/results/TOTTTEND'],nIter));

    %%% Tendency
    tt_tend = squeeze(TOTTTEND)/86400; %%% temperature tendency, degC/s
    tt_tend(tt_tend==0)=NaN;
    
    wN = ww*N2_const*cos_slope/gravity/tAlpha;
    uN = uu*N2_const*sin_slope/gravity/tAlpha;
    
    dtdx = (tt([Nx 1:Nx-1],:) - tt)/dx; %%% on u-grid
    udtdx = uu.*dtdx;
    
    dtdz = zeros(Nx,Nr);
    dtdz(:,2:Nr) = -diff(tt,1,2)./(-diff(zz)); %%% on w-grid
    wdtdz = ww.*dtdz;

    %%% Advection
    dutdx = zeros(Nx,Nr);
    dwtdz = zeros(Nx,Nr);
    
    dutdx = ADVx_TH([2:Nx 1],:) - ADVx_TH(1:Nx,:);
    dwtdz(:,2:Nr) = ADVr_TH(:,2:Nr) - ADVr_TH(:,1:Nr-1);
    tt_adv = dutdx + dwtdz;

    dutdx(dutdx==0)=NaN;
    dutdx = dutdx ./ Cell_Vol;
    dutdx = squeeze(dutdx);
    
    dwtdz(dwtdz==0)=NaN;
    dwtdz = dwtdz ./ Cell_Vol;
    dwtdz = squeeze(dwtdz);
    
    tt_adv(tt_adv==0)=NaN;
    tt_adv = tt_adv ./ Cell_Vol;
    tt_adv = squeeze(tt_adv);
  
    tt = tt + tt_background;
    rho = rhoConst.*(1-(tt-tRef)*tAlpha);
    N2 = NaN*zeros(Nx,Nr);
    N2(:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));


    xidx=1:Nx;
    dutdx_int(o)=sum(dutdx(xidx,Nr-Nbot:Nr),'all','omitnan');
    dwtdz_int(o)=sum(dwtdz(xidx,Nr-Nbot:Nr),'all','omitnan');
    tt_adv_int(o)=sum(tt_adv(xidx,Nr-Nbot:Nr),'all','omitnan');
    tt_tend_int(o)=sum(tt_tend(xidx,Nr-Nbot:Nr),'all','omitnan');
    uN_int(o)=sum(uN(xidx,Nr-Nbot:Nr),'all','omitnan');
    wN_int(o)=sum(wN(xidx,Nr-Nbot:Nr),'all','omitnan');

    % %%% Plot the 4 terms
    % figure(1)
    % clf;set(gcf,'color','w','Position',[101 63 1051 788]);
    % subplot(3,2,1)
    % % pcolor(xx/1000,-zz,udtdx');
    % pcolor(xx/1000,-zz,dutdx');
    % colormap(redblue);
    % hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % % title(['$u\,\partial_x\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % title(['$\partial_x(u\theta)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim(CLIM);ylim(YLIM);xlim(XLIM);
    % 
    % subplot(3,2,2)
    % % pcolor(xx/1000,-zz,wdtdz');
    % pcolor(xx/1000,-zz,dwtdz');
    % hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % % title(['$w\,\partial_z\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % title(['$\partial_z(w\theta)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim(CLIM);ylim(YLIM);xlim(XLIM);
    % 
    % subplot(3,2,3)
    % pcolor(xx/1000,-zz,uN');colormap(redblue);
    % hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$u N^2\sin\theta/(g\alpha)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim(CLIM);ylim(YLIM);xlim(XLIM);
    % 
    % subplot(3,2,4)
    % pcolor(xx/1000,-zz,wN');colormap(redblue);
    % hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$w N^2\cos\theta/(g\alpha)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim(CLIM);ylim(YLIM);xlim(XLIM);
    % 
    % % subplot(3,2,5) %%% Temperature tendency
    % % subplot(3,2,6) %%% Diffusion
    % 
    % subplot(3,2,5)
    % pcolor(xx/1000,-zz,N2')
    % shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$N^2\ (s^{-2})$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim([-2 2]/1e6)
    % ylim(YLIM)
    % 
    % subplot(3,2,6)
    % pcolor(xx/1000,-zz,tt_tend')
    % shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$\partial_t \theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim([-2 2]/1e5)
    % ylim(YLIM)

end

%%

    dutdx_cumint=cumsum(dutdx_int);
    dwtdz_cumint=cumsum(dwtdz_int);
    tt_adv_cumint=cumsum(tt_adv_int);
    tt_tend_cumint=cumsum(tt_tend_int);
    uN_cumint=cumsum(uN_int);
    wN_cumint=cumsum(wN_int);

    % TOT = uN_cumint+wN_cumint+dutdx_cumint+dwtdz_cumint+tt_tend_cumint;
    % TOT = wN_cumint+dutdx_cumint+dwtdz_cumint;
figure(2)
clf;set(gcf,'color','w','Position', [61 280 879 489]);
plot(uN_cumint,'LineWidth',2);
hold on;
plot(1e4*wN_cumint,'LineWidth',2);
plot(1e16*dutdx_cumint,'LineWidth',2);
plot(100*dwtdz_cumint,'LineWidth',2);
% plot(TOT,'LineWidth',2,'Color',gray);
grid on;grid minor;
xlabel('Hours')
set(gca,'Fontsize',fontsize);
l1 = legend('$u N^2\sin\theta/(g\alpha)$',...
    '$10^4\times w N^2\cos\theta/(g\alpha)$',...
    '$10^{16}\times\partial_x(u\theta)$',...
    '$100\times \partial_z(w\theta)$',...
    ...% 'Total advection',...
    'Fontsize',fontsize+3,'interpreter','latex',...
    'Position', [0.1481 0.1213 0.2057 0.2174]);

