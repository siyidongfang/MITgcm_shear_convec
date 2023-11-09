
%%% N^2 equation: d(db/dz)/dt

%%% Compare three terms in the temperature equation:
%%% u*dT/dx, v*dT/dy(==0), w*dT/dz, w*N^2*cos(theta)/g/alpha, u*N^2*sin(theta)/g/alpha

%%% Use ADVr_TH and ADVx_TH, plot phase average over 10 tidal cycles

%%%% TO DO: plot temperature tendency, diffusion, and the total??


clear;close all;
ne =4;
load_all

cons_t2b = gravity*tAlpha; %%% The constant ratio between buoyancy and temperature

o1 = 1;
o2 = 120;
% No = o2-o1+1;
No = 24

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

CLIM = [-1 1]/5e4/20;
YLIM = [1100 1500];XLIM = [-Lx/2/1000 Lx/2/1000];

Hz = sum(delR);
N2const = (1e-3)^2;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);

for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end

dx = delX(1);

%%% Predicted tidal velocity, include the velocity shear
utide = zeros(No,Nx,Nr);
u_shear = 9e-4;
h_shear = 250;
amp_tide = u_shear * h_shear;
pi_const = 3.141592653589793;
omega_tide = 2*pi_const/43200;
shearProfile = ones(1,Nr);
for k=1:Nr
    if((zz(k)-zz(Nr))<h_shear)
        shearProfile(k) = (zz(k)-zz(Nr))/h_shear;
    end
end

for o=1:No
    nIter = dumpIters(o);
    time_s = nIter.*deltaT; %%% time in seconds
    for i=1:Nx
        utide(o,i,:) = amp_tide*cos(omega_tide*time_s)*shearProfile;
    end
end

%%% Check the utide
% aaa = squeeze(utide(:,20,:));
% pcolor(aaa');axis ij;shading flat;colormap(redblue);colorbar;clim([-0.3 0.3]);

dz_tt_tend = zeros(Nx,Nr); 
dz_wN = zeros(Nx,Nr);
dz_uN = zeros(Nx,Nr);
dz_uN_tide = zeros(Nx,Nr);
dz_udtdx = zeros(Nx,Nr);
dz_udtdx_tide = zeros(Nx,Nr);
dz_wdtdz = zeros(Nx,Nr);
dz_dutdx = zeros(Nx,Nr);
dz_dwtdz = zeros(Nx,Nr);
dz_tt_adv = zeros(Nx,Nr);

%%% Integrate over the bottom 250m (?)
Nbot = round(250/delR(end));
dz_udtdx_int = size(1,No);
dz_udtdx_tide_int = size(1,No);
dz_wdtdz_int = size(1,No);
dz_dutdx_int = size(1,No);
dz_dwtdz_int = size(1,No);
dz_tt_adv_int = size(1,No);
dz_tt_tend_int = size(1,No);
dz_uN_int = size(1,No);
dz_uN_tide_int = size(1,No);
dz_wN_int = size(1,No);

for o=1:No
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
    uu_tide = squeeze(utide(o,:,:)); 

    %%% Tendency
    tt_tend = squeeze(TOTTTEND)/86400; %%% temperature tendency, degC/s
    tt_tend(tt_tend==0)=NaN;
   
    wN = ww*N2_const*cos_slope/gravity/tAlpha;
    uN = uu*N2_const*sin_slope/gravity/tAlpha;
    uN_tide = uu_tide*N2_const*sin_slope/gravity/tAlpha;

    dtdx = (tt([Nx 1:Nx-1],:) - tt)/dx; %%% on u-grid
    udtdx = uu.*dtdx;
    udtdx_tgrid = 0.5*(udtdx([2:Nx 1],:)+udtdx);
    udtdx_tide = uu_tide.*dtdx;
    udtdx_tide_tgrid = 0.5*(udtdx_tide([2:Nx 1],:)+udtdx_tide);
    
    dtdz = zeros(Nx,Nr);
    dtdz(:,2:Nr) = -diff(tt,1,2)./(-diff(zz)); %%% on w-grid
    wdtdz = ww.*dtdz;
    wdtdz_center = zeros(Nx,Nr);
    wdtdz_center(:,1:Nr-1) = 0.5*(wdtdz(:,2:Nr)+wdtdz(:,1:Nr-1));

    %%% Advection
    dwtdz = zeros(Nx,Nr);
    
    dutdx = ADVx_TH([2:Nx 1],:) - ADVx_TH(1:Nx,:); %%% t-grid, cell center
    dwtdz(:,2:Nr) = ADVr_TH(:,2:Nr) - ADVr_TH(:,1:Nr-1); %%% t-grid, cell center
    tt_adv = dutdx + dwtdz; %%% t-grid, cell center

    dutdx(dutdx==0)=NaN;
    dutdx = dutdx ./ Cell_Vol;
    dutdx = squeeze(dutdx);
    
    dwtdz(dwtdz==0)=NaN;
    dwtdz = dwtdz ./ Cell_Vol;
    dwtdz = squeeze(dwtdz);
    
    tt_adv(tt_adv==0)=NaN;
    tt_adv = tt_adv ./ Cell_Vol;
    tt_adv = squeeze(tt_adv);
  
    %%%%% Calculate the vertical derivatives
    dz_tt_tend(:,2:Nr)=-diff(tt_tend,1,2)./delR(2:Nr); %%% t-gird, upper level
    dz_wN(:,1:Nr-1)=-diff(wN,1,2)./delR(2:Nr); %%% t-grid, cell center
    dz_uN(:,2:Nr)=-diff(uN,1,2)./delR(2:Nr);   %%% u-grid, upper level
    dz_uN_tide(:,2:Nr)=-diff(uN_tide,1,2)./delR(2:Nr);   %%% u-grid, upper level
    dz_udtdx(:,2:Nr)=-diff(udtdx_tgrid,1,2)./delR(2:Nr); %%% t-grid, upper level
    dz_udtdx_tide(:,2:Nr)=-diff(udtdx_tide_tgrid,1,2)./delR(2:Nr); %%% t-grid, upper level
    dz_wdtdz(:,2:Nr)=-diff(wdtdz_center,1,2)./delR(2:Nr); %%% t-grid, upper level
    dz_dutdx(:,2:Nr)=-diff(dutdx,1,2)./delR(2:Nr); %%% t-grid, upper level
    dz_dwtdz(:,2:Nr)=-diff(dwtdz,1,2)./delR(2:Nr); %%% t-grid, upper level
    dz_tt_adv(:,2:Nr)=-diff(tt_adv,1,2)./delR(2:Nr); %%% t-grid, upper level


    tt = tt + tt_background;
    rho = rhoConst.*(1-(tt-tRef)*tAlpha);
    N2 = NaN*zeros(Nx,Nr);
    N2(:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));

    xidx=100;
    zidx = Nr-10; %%% zidx = Nr-Nbot:Nr;
    dz_udtdx_int(o)=sum(dz_udtdx(xidx,zidx),'all','omitnan');
    dz_wdtdz_int(o)=sum(dz_wdtdz(xidx,zidx),'all','omitnan');
    dz_udtdx_tide_int(o)=sum(dz_udtdx_tide(xidx,zidx),'all','omitnan');
    dz_dutdx_int(o)=sum(dz_dutdx(xidx,zidx),'all','omitnan');
    dz_dwtdz_int(o)=sum(dz_dwtdz(xidx,zidx),'all','omitnan');
    dz_tt_adv_int(o)=sum(dz_tt_adv(xidx,zidx),'all','omitnan');
    dz_tt_tend_int(o)=sum(dz_tt_tend(xidx,zidx),'all','omitnan');
    dz_uN_int(o)=sum(dz_uN(xidx,zidx),'all','omitnan');
    dz_uN_tide_int(o)=sum(dz_uN_tide(xidx,zidx),'all','omitnan');
    dz_wN_int(o)=sum(dz_wN(xidx,zidx),'all','omitnan');

    %%% Plot the 4 terms
    figure(1)
    clf;set(gcf,'color','w','Position',[101 63 1400 1000]);
    subplot(4,3,1)
    % pcolor(xx/1000,-zz,dz_udtdx');
    pcolor(xx/1000,-zz,dz_dutdx');
    colormap(redblue);
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$u\,\partial_x\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    title(['$\partial_x(u\theta)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);
    ylim(YLIM);xlim(XLIM);

    subplot(4,3,2)
    % pcolor(xx/1000,-zz,dz_wdtdz');
    pcolor(xx/1000,-zz,dz_dwtdz');
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$w\,\partial_z\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    title(['$\partial_z(w\theta)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);
    ylim(YLIM);xlim(XLIM);

    subplot(4,3,3)
    pcolor(xx/1000,-zz,dz_uN');colormap(redblue);
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$u N^2\sin\theta/(g\alpha)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);ylim(YLIM);xlim(XLIM);

    subplot(4,3,4)
    pcolor(xx/1000,-zz,dz_wN');colormap(redblue);
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$w N^2\cos\theta/(g\alpha)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);ylim(YLIM);xlim(XLIM);

    subplot(4,3,5)
    pcolor(xx/1000,-zz,N2')
    shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$N^2\ (s^{-2})$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-2 2]/1e6)
    ylim(YLIM)

    subplot(4,3,6)
    pcolor(xx/1000,-zz,dz_tt_tend')
    shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$\partial_t \theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM)
    ylim(YLIM)

    subplot(4,3,7)
    pcolor(xx/1000,-zz,dz_udtdx_tide');
    colormap(redblue);
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$u\,\partial_x\theta$ tide, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);
    ylim(YLIM);xlim(XLIM);

    subplot(4,3,8)
    pcolor(xx/1000,-zz,dz_uN_tide');colormap(redblue);
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$u N^2\sin\theta/(g\alpha)$ tide, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);ylim(YLIM);xlim(XLIM);

    subplot(4,3,9)
    pcolor(xx/1000,-zz,dz_udtdx');
    colormap(redblue);
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$u\,\partial_x\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);
    ylim(YLIM);xlim(XLIM);

    subplot(4,3,10)
    pcolor(xx/1000,-zz,dz_wdtdz');
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$w\,\partial_z\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
     clim(CLIM);
    ylim(YLIM);xlim(XLIM);


    subplot(4,3,11)
    pcolor(xx/1000,-zz,dz_udtdx'+dz_wdtdz');
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$u\,\partial_x\theta$+$w\,\partial_z\theta$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);
    ylim(YLIM);xlim(XLIM);

    subplot(4,3,12)
    pcolor(xx/1000,-zz,dz_dutdx'+dz_dwtdz');
    hold on;shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$\partial_x(u\theta)$+$\partial_z(w\theta)$, t = ' num2str(h_tides,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim(CLIM);
    ylim(YLIM);xlim(XLIM);



end

%%

    dz_dutdx_cumint=cumsum(dz_dutdx_int);
    dz_dwtdz_cumint=cumsum(dz_dwtdz_int);
    dz_tt_adv_cumint=cumsum(dz_tt_adv_int);
    dz_tt_tend_cumint=cumsum(dz_tt_tend_int);
    dz_uN_cumint=cumsum(dz_uN_int);
    dz_wN_cumint=cumsum(dz_wN_int);

    % TOT = uN_cumint+wN_cumint+dutdx_cumint+dwtdz_cumint+tt_tend_cumint;
    % TOT = wN_cumint+dutdx_cumint+dwtdz_cumint;
figure(2)
clf;set(gcf,'color','w','Position', [61 280 879 489]);
plot(dz_uN_cumint,'LineWidth',2);
hold on;
plot(dz_wN_cumint,'LineWidth',2);
plot(dz_dutdx_cumint,'LineWidth',2);
plot(dz_dwtdz_cumint,'LineWidth',2);
plot(dz_tt_tend_cumint,'LineWidth',2);

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

