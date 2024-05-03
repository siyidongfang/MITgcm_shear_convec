%%% 
%%% calc_GrowthRate_MITgcm_RMSE.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations


clear;
% close all

for  ne = 1
load_all

% Ntide = 20;
% tidx = 1:Ntide*12;
No = nDumps;
% No =  round(172286/514)-20;
tidx = 1:No;
Nt = length(tidx);
Hshear = 250;
dz = delR(end);
Nshear = round(Hshear/dz);
zidx = Nr-Nshear:Nr-1;
% zidx = Nr-Nshear+20:Nr-1-20;
% zidx = 1:Nr;
Nshear = length(zidx);
hab_shear = zz(zidx)-min(zz);

div_tt = zeros(Nt,Nshear);
div_uu = zeros(Nt,Nshear);
div_vv = zeros(Nt,Nshear);
div_ww = zeros(Nt,Nshear);

% for o=1:12
%     nIter = dumpIters(o);
%     tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
%     uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
%     vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
%     ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
%     tt_shear = tt(:,zidx);
%     uu_shear = uu(:,zidx);
%     vv_shear = vv(:,zidx);
%     ww_shear = ww(:,zidx);
%     mean_tt_shear12(o,:) = mean(tt(:,zidx));
%     mean_uu_shear12(o,:) = mean(uu(:,zidx));
%     mean_vv_shear12(o,:) = mean(vv(:,zidx));
%     mean_ww_shear12(o,:) = mean(ww(:,zidx));
% end

for o = tidx
    nIter = dumpIters(o);
    time_h(o) = nIter.*deltaT./3600;
    time_tidal(o) = time_h(o)/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    tt_shear = tt(:,zidx);
    uu_shear = uu(:,zidx);
    vv_shear = vv(:,zidx);
    ww_shear = ww(:,zidx);
    
    mean_tt_shear = mean(tt(:,zidx));
    mean_uu_shear = mean(uu(:,zidx));
    mean_vv_shear = mean(vv(:,zidx));
    mean_ww_shear = mean(ww(:,zidx));

    % Nrem = rem(o,12);
    % if(Nrem==0)
    %     Nrem = 12;
    % end
    % mean_tt_shear = mean_tt_shear12(Nrem,:);
    % mean_uu_shear = mean_uu_shear12(Nrem,:);
    % mean_vv_shear = mean_vv_shear12(Nrem,:);
    % mean_ww_shear = mean_ww_shear12(Nrem,:);

    mean_tt_shear_2d = repmat(mean_tt_shear,[Nx 1]);
    mean_uu_shear_2d = repmat(mean_uu_shear,[Nx 1]);
    mean_vv_shear_2d = repmat(mean_vv_shear,[Nx 1]);
    mean_ww_shear_2d = repmat(mean_ww_shear,[Nx 1]);

    div_tt(o,:) = rmse(tt_shear,mean_tt_shear_2d,1);
    div_uu(o,:) = rmse(uu_shear,mean_uu_shear_2d,1);
    div_vv(o,:) = rmse(vv_shear,mean_vv_shear_2d,1);
    div_ww(o,:) = rmse(ww_shear,mean_ww_shear_2d,1);

end

%%

% CLIM = [0 1]*1e-18;
% 
% figure(1)
% clf;set(gcf,'Color','w','Position', [75 224 1362 647])
% subplot(2,2,1)
% pcolor(time_h,hab_shear,div_tt');shading flat;colorbar
% ylabel('HAB (m)');xlabel('Time (hours)');
% title('RMSE of potential temperature (degC)')
% set(gca,'Fontsize',fontsize)
% clim(CLIM)
% subplot(2,2,2)
% pcolor(time_h,hab_shear,div_uu');shading flat;colorbar
% ylabel('HAB (m)');xlabel('Time (hours)');
% set(gca,'Fontsize',fontsize)
% title('RMSE of cross-isobath velocity u RMSE (m/s)')
% clim(CLIM)
% subplot(2,2,3)
% pcolor(time_h,hab_shear,div_vv');shading flat;colorbar
% ylabel('HAB (m)');xlabel('Time (hours)');
% set(gca,'Fontsize',fontsize)
% title('RMSE of along-isobath velocity v (m/s)')
% clim(CLIM)
% subplot(2,2,4)
% pcolor(time_h,hab_shear,div_ww');shading flat;colorbar
% ylabel('HAB (m)');xlabel('Time (hours)');
% set(gca,'Fontsize',fontsize)
% title('RMSE of vertical velocity w (m/s)')
% colormap(WhiteBlueGreenYellowRed(0))
% clim(CLIM)

% % % % % % print('-dpng','-r150',[expname '_rmse.png']);

%%
div_tt2_zavg = mean(div_tt.^2,2);
div_uu2_zavg = mean(div_uu.^2,2);
div_vv2_zavg = mean(div_vv.^2,2);
div_ww2_zavg = mean(div_ww.^2,2);

% div_tt_norm = div_tt_zavg/div_tt_zavg(1);
% div_uu_norm = div_uu_zavg/div_uu_zavg(1);
% div_vv_norm = div_vv_zavg/div_vv_zavg(1);
% div_ww_norm = div_ww_zavg/div_ww_zavg(1);


Pr = 2;
ke = div_uu2_zavg/2+div_vv2_zavg/2+div_ww2_zavg/2;
pe = Pr*div_tt2_zavg/2;
% pe(1:12)=NaN;
% ke(1:12)=NaN;
energy =ke+pe;


% %%% Calculate the growth rate
% fit_span = 12*1+1:6*12;
% if(max(energy)<=1e-5)
%     fit_span = 12*20+1:40*12;
% end
% xxplot = time_h;
% yyplot = log(energy)/2;
% [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% grow(ne) = pp(1)
% [y_fit,delta_fit] = polyval(pp,xxplot,S);

figure()
clf;set(gcf,'Color','w','Position',[211 289 852 394])
plot(time_h/12,log(pe)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(ke)/2,'LineWidth',2);
plot(time_h/12,log(energy)/2,'LineWidth',2);
% plot(xxplot(fit_span)/12, y_fit(fit_span),'--','LineWidth',2);
set(gca,'Fontsize',fontsize)
xlabel('Time (tidal cycles)')
% xlabel('Time (days)')
title('Normalized temperature RMSE')
title('RMSE of T and u averaged over the bottom shear layer')
% title('Vertically averaged RMSE of T and u')
% ylabel('(degC)')
legend('RMSE of PE','RMSE of KE','RMSE of PE+KE','Linear fit')
grid on;grid minor;
hold on;
% ylim([1e-9 1e-1])
% ylim([min(min([div_tt_zavg div_uu_zavg])) max(max([div_tt_zavg div_uu_zavg]))])

% plot(div_uu_norm);
% plot(div_vv_norm);
% plot(div_ww_norm);
print('-dpng','-r150',[expdir expname '_rmse.png']);


% filename = [expdir expname '/RMSE_mean.mat'];

% save(filename,'time_h','hab_shear','div_tt_zavg','div_uu_zavg','div_vv_zavg','div_ww_zavg',...
%     'div_tt_norm','div_vv_norm','div_ww_norm','div_uu_norm',...
%     'div_tt','div_uu','div_vv','div_ww')

% save(filename,'time_h','hab_shear','div_tt_zavg','div_uu_zavg','div_ww_zavg',...
%     'div_tt_norm','div_ww_norm','div_uu_norm',...
%     'div_tt','div_uu','div_ww')


end

