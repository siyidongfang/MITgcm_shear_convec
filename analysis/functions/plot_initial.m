%%%
%%% plot_initial.m
%%%
%%% Plot the initial condition of the simulation


    % For only one experiments
    clear;close all;

    %%% Add path
    addpath /Users/csi/MITgcm_BLT/analysis/;
    addpath /Users/csi/MITgcm_BLT/analysis/functions;
    addpath /Users/csi/MITgcm_BLT/analysis/colormaps;
    addpath /Users/csi/MITgcm_BLT/analysis/colormaps/cmocean/;
    addpath /Users/csi/MITgcm_BLT/analysis/colormaps/customcolormap/;

    load_colors;
    list_exps;
    ne=1;
    expname = EXPNAME{ne}
    loadexp;
    fontsize =18;

    nIter = 0;
    t0 = squeeze(rdmds([exppath,'/results/T'],nIter));
    s0 = squeeze(rdmds([exppath,'/results/S'],nIter));
    u0 = squeeze(rdmds([exppath,'/results/U'],nIter));
    v0 = squeeze(rdmds([exppath,'/results/V'],nIter));
    eta0 = rdmds([exppath,'/results/Eta'],nIter);

    t0(t0==0)=NaN;
    s0(s0==0)=NaN;

    figure(1)
    clf;set(gcf,'color','w','Position',[168 224 1333 326]);
    subplot(1,2,1)
    pcolor(yy/1000,-zz/1000,t0')
    shading flat;colorbar;axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('y (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex')
    title('Initial temperature $(^\circ \mathrm{C})$','Fontsize',fontsize+3,'interpreter','latex')
    colormap(WhiteBlueGreenYellowRed(0));
    subplot(1,2,2)
    pcolor(yy/1000,-zz/1000,s0')
    shading flat;colorbar;axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('y (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex')
    title('Initial salinity (psu)','Fontsize',fontsize+3,'interpreter','latex')
 
