
    clear;close all;
    ne=1;
    load_all;

    t0 = squeeze(rdmds([exppath,'/results/T'],0));

    dumpFreq = abs(diag_frequency(1)); 
    nDumps = floor(nTimeSteps*deltaT/dumpFreq);
    dumpIters = round((1:nDumps)*dumpFreq/deltaT);
    dumpIters = dumpIters(dumpIters > nIter0);

    o2 = 69;
    o1 = 20;

    YLIM = [0 1500];XLIM = [-Lx/2/1000 Lx/2/1000];

    [ZZ,XX] = meshgrid(zz,xx);
    
    Hz = sum(delR);
    N2const = (1e-3)^2;
    tNorth = N2const *(zz+Hz) /9.81/2e-4;
    tt_background = ones(Nx,Nr);

    for k=1:Nr
        tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
    end


for o=o1:o2

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;

    tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
    tt(tt==0)=NaN;
    tt = tt + tt_background;

    uu = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
    uu(uu==0)=NaN;
    ww = squeeze(rdmds([exppath,'/results/WVEL_inst'],nIter));
    ww(ww==0)=NaN;

    eta = rdmds([exppath,'/results/ETAN_inst'],nIter);
    eta(eta==0)=NaN;

    rho = rhoConst.*(1-(tt-tRef)*tAlpha);
    N2 = NaN*zeros(Nx,Nr);
    N2(:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));

    S = NaN*zeros(Nx,Nr);
    S(:,1:Nr-1) = (uu(:,1:end-1)-uu(:,2:end))./(zz(1:end-1)-zz(2:end));

    Ri = NaN*zeros(Nx,Nr);
    Ri = N2./(S.^2);


    %%
    figure(1)
    set(gcf,'Position',  [-407 32 1247 942])
    clf;set(gcf,'color','w');
    subplot(3,2,1)
    pcolor(xx/1000,-zz,tt');hold on;
    % contour(XX/1000,-ZZ,tt,[0:0.01:1],'k')
    shading flat;colorbar;axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['Temperature $\theta (^\circ \mathrm{C})$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-0.1 0.8])
    % clim([-0.1 0.1])
    ylim(YLIM);xlim(XLIM);

    subplot(3,2,2)
    pcolor(xx/1000,-zz,uu')
    shading flat;colorbar;colormap(jet);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['u (m/s), t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-0.3 0.3])
    ylim(YLIM);xlim(XLIM);

    subplot(3,2,3)
    % pcolor(xx/1000,-zz,real(log10(N2))')
    pcolor(xx/1000,-zz,N2')
    shading flat;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$\log(N^2)\ (s^{-2})$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    title(['$N^2\ (s^{-2})$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-2 2]/1e6)
    ylim(YLIM);xlim(XLIM);

    subplot(3,2,4)
    pcolor(xx/1000,-zz,S')
    shading flat;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$du/dz\ (s^{-1})$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-1 1]/1e3)
    ylim(YLIM);xlim(XLIM);

    % subplot(3,2,5)
    % plot(xx/1000,eta,'LineWidth',2);
    % grid on;grid minor;
    % set(gca,'Fontsize',fontsize);
    % xlabel('x (km)','interpreter','latex');
    % title(['Sea Surface Height (m), t = ' num2str(time_h,'%.1f') ' h'],'interpreter','latex');
    % ylim([-0.1 0.1])

    subplot(3,2,5)
    pcolor(xx/1000,-zz,ww')
    shading flat;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['w (m/s), t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-2 2]/100)
    ylim(YLIM);xlim(XLIM);

    subplot(3,2,6)
    pcolor(xx/1000,-zz,1./Ri')
    shading flat;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['Inverse Richardson Number, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-10 10])
    ylim(YLIM);xlim(XLIM);

    % set(gcf, 'InvertHardcopy', 'off')
    % print('-djpeg','-r150',[figdir '/uvtn_' num2str(o) '.jpeg']);
end




