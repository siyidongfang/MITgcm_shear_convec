
clear;close all;
ne =5;
load_all

% Lx
%%

HAB = 100;
n_select = Nr - round(HAB/3);
xx = xx-xx(1)+delX(1)/2;

lambda = diff(xx([1 Nx]))./(0:1:Nx/2-1);
kx = 2*pi./lambda;


 for  o =  110
% for  o = 121-6:121+6
% for o=nDumps-12:nDumps-1
% for o=246-240+1:12:nDumps
    o
    nIter = dumpIters(o);
    tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter)) ;
    tt_selected = tt(:,n_select)';
    detrended_tt= detrend(tt_selected, 1); %%% remove the linear trend
    tt_fft = fft(detrended_tt)/Nx;
    
    fontsize = 16;
    figure()
    set(gcf,'Color','w','Position', [128 250 1153 421])
    subplot(1,2,1)
    loglog(lambda,abs(tt_fft(1:Nx/2)).^2,'LineWidth',2);
    xlabel('Wavelength, \lambda (m)');
    ylabel('Spectral power (^oC^2)');
    title('Potential Temperature Spectrum: HAB = '+string(HAB) +' m')
    grid on;grid minor;
    set(gca,'Fontsize',fontsize)

    subplot(1,2,2)
    loglog(kx, abs(tt_fft(1:Nx/2)).^2,'LineWidth',2);
    xlabel('Wavenumber, k_x (rad/m)');
    ylabel('Spectral power (^oC^2)');
    title('Potential Temperature Spectrum: HAB = '+string(HAB) +' m')
    grid on;grid minor;
    set(gca,'Fontsize',fontsize)


end

