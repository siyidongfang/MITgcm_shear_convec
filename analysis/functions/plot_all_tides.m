

for h_tides = 1:12
  load([prodir expname '_tavg_5days_h' num2str(h_tides) '.mat'])

      tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));


    tt(tt==0)=NaN;
    uu(uu==0)=NaN;
    vv(vv==0)=NaN;

    rhoConst = 999.8;

    rho = rhoConst.*(1-(tt-tRef)*tAlpha);
    N2 = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));
    pp_mid = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined

    figure(3)
    set(gcf,'Position',  [-407 32 1247 942])
    clf;set(gcf,'color','w');
    subplot(1,2,1)
    pcolor(yy/1000,-zz,tt')
    ylim([0 2500]);
    % xlim([spongeThickness*1.5 Ly-spongeThickness*1.5]/1000);
    shading flat;colorbar;axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('y (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex')
    title(['$\theta\ (^\circ \mathrm{C})$, t = ' num2str(time_h,'%.1f') 'h'],'Fontsize',fontsize+3,'interpreter','latex')

    subplot(3,2,4)
    pcolor(yy/1000,-zz,vv')
    ylim([0 2500]);
    shading flat;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('y (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['v (m/s), t = ' num2str(time_h,'%.1f') 'h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([-0.3 0.3])

    subplot(3,2,5)
    pcolor(yy/1000,pp_mid,N2')
    ylim([0 2500]);
    shading flat;colorbar;colormap(redblue);axis ij;set(gca,'Fontsize',fontsize);set(gca,'color',gray);
    xlabel('y (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
    title(['$N^2\ (s^{-2})$, t = ' num2str(time_h,'%.1f') 'h'],'Fontsize',fontsize+3,'interpreter','latex')
    clim([0 2]/1e6)

    set(gcf, 'InvertHardcopy', 'off')

end

