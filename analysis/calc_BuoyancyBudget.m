%%%
%%% calc_BuoyancyBudget.m
%%%
%%% Calculate the buoyancy budget (i.e., temperature budget) for
%%% the 2D BLT simulations

clear;close all;
ne = 1;load_all;

%%% Load data
load([prodir expname '_tavg_5days.mat'])

%%% Grid spacing matrices
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

DZC = zeros(Nx,Ny,Nr);
DZC(:,:,2:Nr) = repmat(reshape(-diff(zz),[1 1 Nr-1]),[Nx Ny 1]);
Cell_Vol = delX(1).*delY(1).*DZC;

%%% Tendency
tt_tend = squeeze(TOTTTEND)/86400; %%% temperature tendency, degC/s
tt_tend(tt_tend==0)=NaN;

%%% Advection
tt_adv = zeros(Nx,Ny,Nr);
dutdx = zeros(Nx,Ny,Nr);
dwtdz = zeros(Nx,Ny,Nr);

j=1;
for k=1:Nr-1
  dutdx(:,j,k+1) = ADVx_TH([2:Nx 1],j,k+1) - ADVx_TH(1:Nx,j,k+1);
  dwtdz(:,j,k+1) = ADVr_TH(:,j,k+1) - ADVr_TH(:,j,k);
  tt_adv(:,j,k+1) = dutdx(:,j,k+1) + dwtdz(:,j,k+1);
end

dutdx(dutdx==0)=NaN;
dutdx = dutdx ./ Cell_Vol;
dutdx = squeeze(dutdx);

dwtdz(dwtdz==0)=NaN;
dwtdz = dwtdz ./ Cell_Vol;
dwtdz = squeeze(dwtdz);

tt_adv(tt_adv==0)=NaN;
tt_adv = tt_adv ./ Cell_Vol;
tt_adv = squeeze(tt_adv);


%%% Horizontal diffusion
tt_hori = zeros(Nx,Ny,Nr);
i=1;
for j=1:Ny-1
    for k=1:Nr-1
      tt_hori(i,j,k+1) = DFyE_TH(i,j,k) - DFyE_TH(i,j+1,k); 
    end
end

tt_hori(tt_hori==0)=NaN;
tt_hori = tt_hori ./ Cell_Vol;
tt_hori = squeeze(tt_hori);


%%% Vertical mixing: implicit part
tt_vertImpl = zeros(Nx,Ny,Nr);
i=1;
for j=1:Ny-1
    for k=1:Nr-1
      tt_vertImpl(i,j,k+1) = DFrI_TH(i,j,k+1) - DFrI_TH(i,j,k);
    end
end
tt_vertImpl(tt_vertImpl==0)=NaN;
tt_vertImpl = tt_vertImpl ./ Cell_Vol;
tt_vertImpl = squeeze(tt_vertImpl);

% tt_vert_difference = tt_vertImpl - tt_vertExp; %%% almost zero

%%% Adams-Bashforth --- always zero
% % % AB_gT   | 15 |       |SMR     MR|degC/s          |Potential Temp. tendency from Adams-Bashforth
% % % gTinAB  | 15 |       |SMR     MR|degC/s          |Potential Temp. tendency going in Adams-Bashforth
% tt_ab = squeeze(AB_gT);
% tt_ab(tt_ab==0)=NaN;

%%% Residual
% tt_residual = tt_adv + tt_vertImpl + tt_ab - tt_tend;
tt_residual = tt_adv + tt_hori + tt_vertImpl - tt_tend;


[ZZ,YY] = meshgrid(zz,yy);
CLIM = [-1 1]/1e7;
figure(1)
clf;
set(gcf,'color','w'); set(gcf,'Position',[136 144 1158 1000]);
subplot(3,2,1)
pcolor(YY/1000,-ZZ/1000,tt_adv); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim([1.1 -min(bathy)/1000+0.05]);
% xlim([0.75 14.25])
title('$\theta$ budget: Advection ($^o$C/s)','interpreter','latex','FontSize',fontsize+4);

subplot(3,2,2)
pcolor(YY/1000,-ZZ/1000,tt_tend); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim([1.1 -min(bathy)/1000+0.05]);
% xlim([0.75 14.25])
title('$\theta$ budget: Tendency ($^o$C/s)','interpreter','latex','FontSize',fontsize+4);

subplot(3,2,3)
pcolor(YY/1000,-ZZ/1000,tt_vertImpl); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/100);
ylim([1.1 -min(bathy)/1000+0.05]);
% xlim([0.75 14.25])
title('$\theta$ budget: Vertical mixing ($^o$C/s)','interpreter','latex','FontSize',fontsize+4);
 
subplot(3,2,4)
pcolor(YY/1000,-ZZ/1000,tt_hori); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/100);
ylim([1.1 -min(bathy)/1000+0.05]);
% xlim([0.75 14.25])
title('$\theta$ budget: Horizontal diffusion ($^o$C/s)','interpreter','latex','FontSize',fontsize+4);

subplot(3,2,5)
pcolor(YY/1000,-ZZ/1000,tt_residual); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim([1.1 -min(bathy)/1000+0.05]);
% xlim([0.75 14.25])
title('$\theta$ budget: Residual ($^o$C/s)','interpreter','latex','FontSize',fontsize+4);



% print('-djpeg','-r150',[figdir 'BuoyancyBudget.jpeg']);




