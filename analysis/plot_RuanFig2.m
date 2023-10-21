%%%
%%% plot_RuanFig2.m
%%%
%%% Replicate Fig. 2 of Xiaozhou's manuscript

clear;close all;
ne = 1;
load_all


No = 360; 
vv_timeseries = zeros(No,Nr);
% N2_timeseries = zeros(No,Nr-1);
N2_timeseries = zeros(No,Nr);
time_tidal = zeros(1,No);
pp_mid = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined

dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

tRef = 0;

for o=1:No

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    time_tidal(o) = time_h/12;

    tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
    % ss = squeeze(rdmds([exppath,'/results/SALT'],nIter));
    vv = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
    % n2 = -gravity/rhoConst*squeeze(rdmds([exppath,'/results/DRHODR'],nIter));

    tt(tt==0)=NaN;
    % ss(ss==0)=NaN;
    vv(vv==0)=NaN;
    % n2(n2==0)=NaN;

    % LC = 0.75*Ly;
    LC = 0.85*Ly;
    [mC,nC] = min(abs(yy-LC));
    ttC = tt(nC,:);
    % ssC = tt(nC,:);
    vv_timeseries(o,:) = vv(nC,:);
    tt_timeseries(o,:) = tt(nC,:);
    rhoC = rhoConst.*(1-(ttC-tRef)*tAlpha);
    N2_timeseries(o,2:Nr) = -gravity/rhoConst.*(rhoC(1:end-1)-rhoC(2:end))./(zz(1:end-1)-zz(2:end));
    % N2_timeseries(o,:) = n2(nC,:);

end


botN = find(isnan(vv_timeseries(1,:)),1)
botZ = zz(botN)

%%
figure(1)
set(gcf,'Position', [139 266 833 623])
clf;set(gcf,'color','w');
subplot(2,1,1)
pcolor(time_tidal,zz-botZ,vv_timeseries');
hold on;
contour(time_tidal,zz-botZ,vv_timeseries',[-0.16:0.02:0.16],'color','k')
% pcolor(time_tidal,-zz,vv_timeseries')
shading interp;colorbar;colormap(redblue);set(gca,'Fontsize',fontsize);set(gca,'color',gray);
xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('u (m/s)','Fontsize',fontsize+4,'interpreter','latex')
clim([-0.2 0.2])
ylim([-30 500])

subplot(2,1,2)
pcolor(time_tidal,zz-botZ,(N2_timeseries)')
hold on;contour(time_tidal,zz-botZ,(N2_timeseries)',[0 0],'Color',cyan,'LineWidth',2);
% pcolor(time_tidal,-pp_mid-botZ,N2_timeseries')
% pcolor(time_tidal,zz-botZ,real(log10(sqrt(N2_timeseries)))')
% pcolor(time_tidal,pp_mid,N2_timeseries')
shading interp;colorbar;
% colormap(cmocean('delta'));
set(gca,'Fontsize',fontsize);set(gca,'color',gray);
xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('$N^2\ (s^{-2})$','Fontsize',fontsize+4,'interpreter','latex')
clim([-3 3]/1e6)
% clim([0 3]/1e6)
 % clim([0.98 1.02]/1e6)
ylim([-30 300])

figure(3)
pcolor(time_tidal,zz-botZ,tt_timeseries')
shading interp;colorbar;
set(gca,'Fontsize',fontsize);set(gca,'color',gray);
xlabel('Tidal cycles','interpreter','latex');ylabel('HAB (m)','interpreter','latex')
title('$t$','Fontsize',fontsize+4,'interpreter','latex')
clim([0.2 0.6]);
colormap(redblue)
ylim([-30 300])
 set(gcf, 'InvertHardcopy', 'off')
