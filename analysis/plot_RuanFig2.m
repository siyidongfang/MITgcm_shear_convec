%%%
%%% plot_RuanFig2.m
%%%
%%% Replicate Fig. 2 of Xiaozhou's manuscript

clear;close all;
ne =1;
load_all

% xx = xx-xx(1);

No = 319; 
uu_timeseries = zeros(No,Nr);
% vv_timeseries = zeros(No,Nr);
N2_timeseries = zeros(No,Nr);
time_tidal = zeros(1,No);
pp_mid = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined

dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

tRef = 0;
m1km = 1000;

    Hz = sum(delR);
    N2const = (1e-3)^2;
    tNorth = N2const *(zz+Hz) /9.81/2e-4;
    tt_background = ones(Nx,Nr);

    for k=1:Nr
        tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
    end


for o=1:No

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    time_tidal(o) = time_h/12;

    tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
    % ss = squeeze(rdmds([exppath,'/results/SALT'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL_inst'],nIter));
    % n2 = -gravity/rhoConst*squeeze(rdmds([exppath,'/results/DRHODR'],nIter));

    tt(tt==0)=NaN;
    tt = tt+tt_background;

    % ss(ss==0)=NaN;
    uu(uu==0)=NaN;
    vv(vv==0)=NaN;
    % n2(n2==0)=NaN;

    % LC = 0.75*Ly;
    % LC = 0.85*Lx;

    LC = 2*m1km;
    [mC,nC] = min(abs(xx-LC));
    ttC = tt(nC,:);
    % ssC = tt(nC,:);
    uu_timeseries(o,:) = uu(nC,:);
    vv_timeseries(o,:) = vv(nC,:);
    tt_timeseries(o,:) = tt(nC,:);
    rhoC = rhoConst.*(1-(ttC-tRef)*tAlpha);
    N2_timeseries(o,2:Nr) = -gravity/rhoConst.*(rhoC(1:end-1)-rhoC(2:end))./(zz(1:end-1)-zz(2:end));
    % N2_timeseries(o,:) = n2(nC,:);

end


botN = find(isnan(uu_timeseries(1,:)),1)
% botZ = zz(botN)
botZ =-1500;

%%

% YLIM = [-30 300];
YLIM = [0 600];
% YLIM = [0 1500];
XLIM = [0 25];

figure(1)
set(gcf,'Position',[56 139 898 762])
clf;set(gcf,'color','w');
subplot(3,1,1)
pcolor(time_tidal,zz-botZ,uu_timeseries');
hold on;
contour(time_tidal,zz-botZ,uu_timeseries',[0.1:0.1:1],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.1],'--','color',darkgray)
% contour(time_tidal,zz-botZ,uu_timeseries',[-0.15:0.03:0.15],'color','k')
shading interp;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('u (m/s)','Fontsize',fontsize+4,'interpreter','latex')
clim([-0.7 0.7])
ylabel('HAB (m)','interpreter','latex')
ylim(YLIM)
xlim(XLIM)

subplot(3,1,2)
pcolor(time_tidal,zz-botZ,tt_timeseries')
hold on;
shading interp;colorbar;
contour(time_tidal,zz-botZ,uu_timeseries',[0.1:0.1:1],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.1],'--','color',darkgray)
set(gca,'Fontsize',fontsize);set(gca,'color',gray);
% xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('$\theta \ (^\circ \mathrm{C})$','Fontsize',fontsize+4,'interpreter','latex')
ylabel('HAB (m)','interpreter','latex')
clim([-0.1 0.8]);
% clim([-0.1 0.3]);
% clim([-0.1 0.1]);
colormap(redblue)
ylim(YLIM)
xlim(XLIM)

subplot(3,1,3)
pcolor(time_tidal,zz-botZ,(N2_timeseries)')
hold on;contour(time_tidal,zz-botZ,(N2_timeseries)',[0 0],'Color',cyan,'LineWidth',3);
% pcolor(time_tidal,-pp_mid-botZ,N2_timeseries')
% pcolor(time_tidal,zz-botZ,real(log10(sqrt(N2_timeseries)))')
% pcolor(time_tidal,pp_mid,N2_timeseries')
contour(time_tidal,zz-botZ,uu_timeseries',[0.1:0.1:1],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1.5)
contour(time_tidal,zz-botZ,uu_timeseries',[-1:0.1:-0.1],'--','color',darkgray)
shading interp;colorbar;
% colormap(cmocean('delta'));
set(gca,'Fontsize',fontsize);set(gca,'color',gray);
xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('$N^2\ (s^{-2})$','Fontsize',fontsize+4,'interpreter','latex')
clim([-3 3]/1e6)
% clim([0 3]/1e6)
 % clim([0.98 1.02]/1e6)
ylim(YLIM)
xlim(XLIM)


figure(2)
pcolor(time_tidal,zz-botZ,vv_timeseries');
hold on;
contour(time_tidal,zz-botZ,vv_timeseries',[-1:0.1:1],'color','k')
shading interp;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('v (m/s)','Fontsize',fontsize+4,'interpreter','latex')
clim([-1 1])
ylim(YLIM)
xlim(XLIM)


 set(gcf, 'InvertHardcopy', 'off')
