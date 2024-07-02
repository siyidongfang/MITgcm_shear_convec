
clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

%--- load the data
addpath ../analysis/
addpath ../analysis/functions/
% expname = 'topo0_H500_s0.0017dz1dx3ln200n-20sm100_kv2e-4';
% expdir = '../exps_hires/';
expname = 'hires_topo4_s0.0006_dz1dx6n-20';
expdir = '/Volumes/MIT/MITgcm_shear_convec/exps_topo4_test/';
loadexp;
rhoConst = 999.8;

dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);
Ntide = 28;
No = Ntide*12;

%--- calculate the mean flow of each tidal cycle
uu_mean= zeros(Ntide,Nr);
ww_mean= zeros(Ntide,Nr);

for nt = 1:Ntide
    nt
    ntidx = nt*12+1:nt*12+12;
    uu = zeros(1,Nr);
    ww = zeros(1,Nr);
    for o=ntidx
        nIter = dumpIters(o);
        time_h = nIter.*deltaT./3600;
        time_tidal(o) = time_h/12;
        uu = uu+mean(squeeze(rdmds([exppath,'/results/UVEL'],nIter)));  
        ww = ww+mean(squeeze(rdmds([exppath,'/results/WVEL'],nIter)));  
    end
    uu_mean(nt,:) = uu/12;
    ww_mean(nt,:) = ww/12;
end

%%
umean = mean(uu_mean(15:end,:));
wmean = mean(ww_mean(15:end,:));

YLIM = [0 500];
fontsize=18;

tt_tide =1:Ntide;
botZ =zz(end);

figure(1)
clf;set(gcf,'Color','w')
subplot(2,2,1)
pcolor(tt_tide,zz-botZ,uu_mean');shading flat;colorbar;
colormap(redblue);clim([-0.1 0.1]/10)
ylabel('HAB (m)','interpreter','latex')
title({'Mean horizontal velocity everaged','over each tidal cycle (m/s)'})
set(gca,'FontSize',fontsize)
xlabel('Number of tidal cycles')
ylim(YLIM)

subplot(2,2,2)
plot(umean,zz-botZ,'LineWidth',2)
hold on;
plot(0.*umean,zz-botZ,'k--')
ylabel('HAB (m)','interpreter','latex')
title({'Mean horizontal velocity averaged','over multiple tidal cycles after convection'})
xlabel('(m/s)')
grid on;
set(gca,'FontSize',fontsize)
ylim(YLIM)


subplot(2,2,3)
pcolor(tt_tide,zz-botZ,ww_mean');shading flat;colorbar;
colormap(redblue);clim([-0.1 0.1]/4e17*1e10)
ylabel('HAB (m)','interpreter','latex')
title({'Mean vertical velocity everaged','over each tidal cycle (m/s)'})
set(gca,'FontSize',fontsize)
xlabel('Number of tidal cycles')
ylim(YLIM)

subplot(2,2,4)
plot(wmean,zz-botZ,'LineWidth',2)
hold on;
plot(0.*wmean,zz-botZ,'k--')
ylabel('HAB (m)','interpreter','latex')
title({'Mean vertical velocity averaged','over multiple tidal cycles after convection'})
xlabel('(m/s)')
grid on;
set(gca,'FontSize',fontsize)
ylim(YLIM)

%--- save the data
