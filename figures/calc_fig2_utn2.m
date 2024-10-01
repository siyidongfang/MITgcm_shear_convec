addpath ../analysis/
addpath ../analysis/functions/
% expname = 'topo0_H500_s0.0017dz1dx3ln200n-20sm100_kv2e-4';
% expdir = '../exps_hires/';
% expname = 'hires_topo4_s0.0013_dz1dx6n-20';
% expdir = '/Volumes/MIT/MITgcm_shear_convec/exps_topo4_test/';

% %%% old
% expname = 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
% expdir = '../exps_topo4/';

%%% new
% expdir = '../exps_kv5e-6/';
% expname='topo4_H500Lx3k_s1.70dz1dx3sm100'
expdir = '../exps_kv1e-5/';
expname='topo4_H500Lx3k_s1.700dz1dx3sm100'
topo = 4;

% expname='flat_H500Lx3k_s1.80dz1dx3sm100';
% topo = 0;

% expdir = '/Volumes/MIT/MITgcm_shear_convec/exps_topo4/';

% expdir = '../exps_topo4_hires/';
% expname = 'topo4_H500Lx3k_s1.6dz1dx3n-20sm100_kv1e-4';

loadexp;

rhoConst = 999.8;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


% xx = xx-xx(1);
% No = nDumps;
% No = 30*12;
No = 28*12;
uu_timeseries= zeros(No,Nr);
% shear_timeseries = zeros(No,Nr);
% vv_timeseries = zeros(No,Nr);
% ww_timeseries = zeros(No,Nr);
tt_timeseries = zeros(No,Nr);
N2_timeseries = zeros(No,Nr);


time_tidal = zeros(1,No);
pp_mid = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined


tRef = 0;
m1km = 1000;

Hz = sum(delR);
N2const = 1e-6;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);

for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end


for o=1:No

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    time_tidal(o) = time_h/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    tt(tt==0)=NaN;

    LC = 0.3*m1km;
    [mC,nC] = min(abs(xx-LC));
    ttC = tt(nC,:)+tt_background(nC,:);
    tt_timeseries(o,:) = tt(nC,:);
    rhoC = rhoConst.*(1-(ttC-tRef)*tAlpha);
    N2_timeseries(o,2:Nr) = -cosd(topo)*gravity/rhoConst.*(rhoC(1:end-1)-rhoC(2:end))./(zz(1:end-1)-zz(2:end));


    ttC1 = tt(nC-1,:)+tt_background(nC,:);
    ttC3 = tt(nC+1,:)+tt_background(nC,:);
    rhoC1 = rhoConst.*(1-(ttC1-tRef)*tAlpha);
    rhoC3 = rhoConst.*(1-(ttC3-tRef)*tAlpha);
    rho_ugrid1 = 0.5*(rhoC1+rhoC);
    rho_ugrid2 = 0.5*(rhoC+rhoC3);
    N2x = -sind(topo)*gravity/rhoConst.*(rho_ugrid2-rho_ugrid1)./delX(1);
    N2x_lower = NaN*zeros(1,Nr);
    N2x_lower(1:Nr-1) = 0.5*(N2x(1:end-1)+N2x(2:end));

    N2_timeseries(o,:) = N2_timeseries(o,:) + N2x_lower;
    

    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    % vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    uu(uu==0)=NaN;
    ww(ww==0)=NaN;
    % vv(vv==0)=NaN;
    shear_timeseries(o,2:Nr) = (uu(nC,1:end-1)-uu(nC,2:end))./delR(2:end);
    uu_timeseries(o,:) = uu(nC,:);
    % vv_timeseries(o,:) = vv(nC,:);
    ww_timeseries(o,:) = ww(nC,:);

end


botN = Nr;
botZ =zz(end);

save('fig2/fig2_new_1.7_kv1e-5.mat','time_tidal','zz','botZ','tt_timeseries','N2_timeseries','uu_timeseries')

% save('figS_gcm_flat/figS_gcm_flat_1.8.mat','time_tidal','zz','botZ','tt_timeseries','N2_timeseries','uu_timeseries')

% %%
% 
% YLIM = [0 Hz];
% XLIM = [0 No/12];
% 
% figure(1)
% set(gcf,'Position',[56 139 898 762])
% clf;set(gcf,'color','w');
% subplot(3,1,2)
% pcolor(time_tidal,zz-botZ,tt_timeseries')
% hold on;
% shading flat;colorbar;
% % contour(time_tidal,zz-botZ,uu_timeseries',[0.05:0.1:1],'color',darkgray)
% % contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
% % contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.05],'--','color',darkgray)
% set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% % xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
% title('$\theta^\prime \ (^\circ \mathrm{C})$','Fontsize',fontsize+4,'interpreter','latex')
% ylabel('HAB (m)','interpreter','latex')
% clim([-0.1 0.1]/5);
% colormap(redblue)
% ylim(YLIM)
% xlim(XLIM)
% 
% subplot(3,1,3)
% pcolor(time_tidal,zz-botZ,(N2_timeseries)')
% hold on;
% contour(time_tidal,zz-botZ,(N2_timeseries)',[0 0],'Color','c','LineWidth',2);
% % contour(time_tidal,zz-botZ,uu_timeseries',[0.05:0.1:1],'color',darkgray)
% % contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
% % contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.05],'--','color',darkgray)
% shading flat;colorbar;
% % colormap(cmocean('delta'));
% set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
% title('$N^2\ (s^{-2})$','Fontsize',fontsize+4,'interpreter','latex')
% clim(([-1 1]+1)/1e6)
% ylim(YLIM)
% xlim(XLIM)
% 
% subplot(3,1,1)
% pcolor(time_tidal,zz-botZ,uu_timeseries');
% hold on;
% contour(time_tidal,zz-botZ,uu_timeseries',[0.05:0.1:1],'color',darkgray)
% contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
% contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.05],'--','color',darkgray)
% shading flat;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% % xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
% title('u (m/s)','Fontsize',fontsize+4,'interpreter','latex')
% clim([-0.6 0.6])
% ylabel('HAB (m)','interpreter','latex')
% ylim(YLIM)
% xlim(XLIM)
% 
% 
% figure(2)
% set(gcf,'Position',[56 139 898 762])
% clf;set(gcf,'color','w');
% subplot(3,1,1)
% pcolor(time_tidal,zz-botZ,ww_timeseries');
% hold on;
% contour(time_tidal,zz-botZ,uu_timeseries',[0.05:0.1:1],'color',darkgray)
% contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
% contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.05],'--','color',darkgray)
% shading flat;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
% title('w (m/s)','Fontsize',fontsize+4,'interpreter','latex')
% clim([-1 1]/100/1e2)
% ylim(YLIM)
% xlim(XLIM)
% 
% 
% subplot(3,1,2)
% pcolor(time_tidal,zz-botZ,vv_timeseries');
% hold on;
% contour(time_tidal,zz-botZ,uu_timeseries',[0.05:0.1:1],'color',darkgray)
% contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
% contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.05],'--','color',darkgray)
% shading interp;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
% title('v (m/s)','Fontsize',fontsize+4,'interpreter','latex')
% clim([-0.3 0.3])
% ylim(YLIM)
% xlim(XLIM)
% 
% 
% subplot(3,1,3)
% pcolor(time_tidal,zz-botZ,shear_timeseries');
% hold on;
% contour(time_tidal,zz-botZ,uu_timeseries',[0.05:0.1:1],'color',darkgray)
% contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
% contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.05],'--','color',darkgray)
% shading interp;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
% title('shear (1/s)','Fontsize',fontsize+4,'interpreter','latex')
% clim([-0.3 0.3]/100)
% ylim(YLIM)
% xlim(XLIM)
% 
% % print('-dpng','-r150',[expdir expname '_fig2.png']);