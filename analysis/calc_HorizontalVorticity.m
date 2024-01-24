%%%
%%% calc_HorizontalVorticity.m
%%%
%%% Calculate zeta = \partial_x w - \partial_z u

clear;close all;
%%% Select an experiment
ne=1;
load_all;

%%% Load data
for o=1:100
    nIter = dumpIters(o);
    uu = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL_inst'],nIter));
    
    DZ = repmat(reshape(squeeze(DRC(2:Nr))',[1 Nr-1]),[Nx 1]);
    
    %%% Calculate vorticity
    dx = delX(1);
    dwdx = diff(ww([Nx 1:Nx],:))/dx;   %%% u-grid, w-level
    dudz = NaN*zeros(Nx,Nr);
    dudz(:,2:Nr) = - diff(uu,1,2)./DZ; %%% u-grid, w-level
    zeta = dwdx - dudz;
    
    figure(1)
    set(gcf,'Position', [117 580 1612 335])
    subplot(1,3,1)
    pcolor(xx/1000,-zz/1000,dwdx')
    shading flat;colorbar;clim([-1 1]/1e5);axis ij;
    subplot(1,3,2)
    pcolor(xx/1000,-zz/1000,-dudz')
    shading flat;colorbar;clim([-1 1]/1e5);axis ij;
    subplot(1,3,3)
    pcolor(xx/1000,-zz/1000,zeta')
    shading flat;colorbar;clim([-1 1]/1e5);axis ij;
    colormap(redblue)

end





