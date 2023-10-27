


 
    clear;
    close all;
    ne=1;
    load_all;
    
    dumpFreq = abs(diag_frequency(1)); 
    nDumps = floor(nTimeSteps*deltaT/dumpFreq);
    dumpIters = round((1:nDumps)*dumpFreq/deltaT);
    dumpIters = dumpIters(dumpIters > nIter0);

    % t0 = squeeze(rdmds([exppath,'/results/T'],0));

    yy=xx;

    o2 = 480;
    o1 = 400;
    % YLIM = [1000 1500]
    YLIM = [0 1500]

    Hz = sum(delR);
    N2const = (1e-3)^2;
    tNorth = N2const *(zz+Hz) /9.81/2e-4;
    tt_background = ones(Nx,Nr);

    for k=1:Nr
        tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
    end


for o=o1:o2

    nIter = dumpIters(o)
    time_h = nIter.*deltaT./3600;
    % nIter = 0

    tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
    tt(tt==0)=NaN;
    tt = tt+tt_background;
    
    vv = squeeze(rdmds([exppath,'/results/VVEL_inst'],nIter));
    vv(vv==0)=NaN;
    uu = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
    uu(uu==0)=NaN;



    eta = rdmds([exppath,'/results/ETAN_inst'],nIter);
    % eta(eta==0)=NaN;


    rho = rhoConst.*(1-(tt-tRef)*tAlpha);
    N2 = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));
    pp_mid = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined

    figure(1)
    set(gcf,'Position',  [-407 32 1247 942])
    clf;set(gcf,'color','w');
    subplot(3,2,1)
    pcolor(yy/1000,-zz,tt')
    % pcolor(yy/1000,-zz,tt'+tt_background')
    % pcolor(yy/1000,-zz,tt'-t0')
    shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['Temperature $\theta (^\circ \mathrm{C})$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % clim([-0.1 0.1])
    % clim([-0.1 0.3])
    clim([-0.1 0.8])
    ylim(YLIM)

    % subplot(3,2,2)
    % pcolor(yy/1000,-zz,ss'-s0')
    % clim([-0.05 0.05])
    % shading interp;colorbar;axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    % xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex')
    % title(['Salinity anomaly $S-S_0$ (psu), t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    % % title('$\theta\ (^\circ \mathrm{C})$','Fontsize',fontsize+3,'interpreter','latex')

    subplot(3,2,3)
    pcolor(yy/1000,-zz,uu')
    shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['u (m/s), t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-0.3 0.3]*3)
    % title('u (m/s)','Fontsize',fontsize+3,'interpreter','latex')
    % clim([-0.3 0.3]/10)
    ylim(YLIM)
    
    subplot(3,2,4)
    pcolor(yy/1000,-zz,vv')
    % ylim([-max(bathy)/1000-0.2 -min(bathy)/1000+0.05]);
    shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['v (m/s), t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-0.3 0.3])
    ylim(YLIM)
    % title('v (m/s)','Fontsize',fontsize+3,'interpreter','latex')

    subplot(3,2,5)
    pcolor(yy/1000,pp_mid,N2')
    % pcolor(yy/1000,-zz/1000,real(log10(N2))')
    % ylim([-max(bathy)/1000-0.2 -min(bathy)/1000+0.05]);
    shading interp;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    % title(['$log(N^2)\ (s^{-2})$, t = ' num2str(time_h,'%.1f') ' h'],'Fontsize',fontsize+3,'interpreter','latex')
    title('$N^2\ (s^{-2})$','Fontsize',fontsize+3,'interpreter','latex')
    % clim([-7 -4])
    clim([-2 2]/1e6)
    ylim(YLIM)

    subplot(3,2,6)
    plot(yy/1000,eta,'LineWidth',2);
    grid on;grid minor;
    set(gca,'Fontsize',fontsize);
    xlabel('x (km)','interpreter','latex');
    title('Sea Surface Height (m)','interpreter','latex');
    ylim([-1 1])
    % ylim([-0.02 0.02])
    

    set(gcf, 'InvertHardcopy', 'off')
    % print('-djpeg','-r150',[figdir '/5dayMean_uvtn.jpeg']);
    % print('-djpeg','-r150',[figdir '/uvtn_' num2str(o) '.jpeg']);
end




