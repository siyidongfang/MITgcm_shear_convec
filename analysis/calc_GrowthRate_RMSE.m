%%% 
%%% calc_GrowthRate_RMSE.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations


clear;
% close all

for  ne = 1
load_all

% Ntide = 20;
% tidx = 1:Ntide*12;
No = nDumps;
No = 35*12;
tidx = 1:No;
Nt = length(tidx);
Hshear = 250;
dz = delR(end);
Nshear = round(Hshear/dz);
% Nshear = 250;
% zidx = Nr-Nshear:Nr;
zidx = 300:500;
Nshear = length(zidx);

div_tt = zeros(Nt,Nshear);
div_uu = zeros(Nt,Nshear);
div_vv = zeros(Nt,Nshear);
div_ww = zeros(Nt,Nshear);

parfor o = tidx
    nIter = dumpIters(o);
    time_h(o) = nIter.*deltaT./3600;
    time_tidal(o) = time_h(o)/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    tt_shear = tt(:,zidx);
    mean_tt_shear = mean(tt(:,zidx));
    mean_tt_shear_2d = repmat(mean_tt_shear,[Nx 1]);
    div_tt(o,:) = rmse(tt_shear,mean_tt_shear_2d,1);

    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    uu_shear = uu(:,zidx);
    vv_shear = vv(:,zidx);
    ww_shear = ww(:,zidx);
    mean_uu_shear = mean(uu(:,zidx));
    mean_vv_shear = mean(vv(:,zidx));
    mean_ww_shear = mean(ww(:,zidx));
    mean_uu_shear_2d = repmat(mean_uu_shear,[Nx 1]);
    mean_vv_shear_2d = repmat(mean_vv_shear,[Nx 1]);
    mean_ww_shear_2d = repmat(mean_ww_shear,[Nx 1]);
    div_uu(o,:) = rmse(uu_shear,mean_uu_shear_2d,1);
    div_vv(o,:) = rmse(vv_shear,mean_vv_shear_2d,1);
    div_ww(o,:) = rmse(ww_shear,mean_ww_shear_2d,1);

end

%%
div_tt2_zint = sum(dz*div_tt.^2,2);
div_uu2_zint = sum(dz*div_uu.^2,2);
div_vv2_zint = sum(dz*div_vv.^2,2);
div_ww2_zint = sum(dz*div_ww.^2,2);

Pr = 1;
ke = (div_uu2_zint/2+div_vv2_zint/2+div_ww2_zint/2);
% ke = 0;
ke = ke/ke(end);
pe = Pr*div_tt2_zint/2;
pe = pe/pe(end);
%%% Calculate the growth rate
% fit_span = 4*12+1:14*12;
% fit_span = 22*12+1:No;
fit_span = 3*12+1:9*12;

% if(max(energy)<=1e-5)
%     % fit_span = 12*10+1:No;
%     fit_span = 12*45+1:No;
% end
xxplot = time_h;
yyplot = log(ke)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(ne) = pp(1)
[y_fit,delta_fit] = polyval(pp,xxplot,S);


% filename = [expdir expname '/RMSE.mat'];
% save(filename,'time_h','xxplot','yyplot','fit_span','pp','y_fit',...
%     'ke','pe','div_uu2_zint','div_vv2_zint','div_ww2_zint','div_tt2_zint','grow')
% % save(filename,'time_h','xxplot','yyplot','fit_span','pp','y_fit',...
% %     'energy','ke','pe','div_tt2_zavg','grow')



h_figure = figure(1);
% set(h_figure,'Visible',false);
clf;
set(gcf,'Color','w','Position',[211 289 852 394])
plot(time_h/12,log(pe)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(ke)/2,'LineWidth',2);
plot(xxplot(fit_span)/12, y_fit(fit_span),'k-.','LineWidth',2);
set(gca,'Fontsize',fontsize+5)
xlabel('Time (tidal cycles)','Interpreter','latex')
% xlabel('Time (days)')
% title('Normalized temperature RMSE')
title('Normalized turbulent energy in the shear layer','Interpreter','latex')
ylabel('$\log(\mathrm{turbulent\ energy})$','Interpreter','latex')
% title('RMSE of T and u averaged over the bottom shear layer')
% title('Vertically averaged RMSE of T and u')
legend('Turbulent potential energy','Turbulent kinetic energy','Linear fit','Position',[0.6972 0.1497 0.1890 0.2069]...
    ,'Interpreter','latex')
grid on;grid minor;
hold on;
% ylim([1e-9 1e-1])
% ylim([min(min([div_tt_zavg div_uu_zavg])) max(max([div_tt_zavg div_uu_zavg]))])

ylim([-40 5])
% print('-dpng','-r150',[expdir expname '_rmse_all.png']);



end

% legend(...
%     's0.0006, NOsmag, h1e-4v2e-4',...
%     's0.0006, NOsmag, h1e-4v5e-4',...
%     's0.0006, NOsmag, h4e-4v2e-4',...
%     's0.0006, smag4C, h1e-4v2e-4',...
%     's0.0006, smag4C, h1e-4v5e-4',...
%     's0.0006, smag4C, h5e-4v2e-4',...
%     's0.0006, 3DSmag1e-2, hv1e-5',...
%     's0.0006, 3Dsmag1e-3, hv1e-5',...
%     's0.0006, 3Dsmag1e-4, hv1e-5')

% legend('smag, hori=1e-4,vert=2e-4','smag, hori=1e-4,vert=5e-4',...
%     'smag, hori=1e-4,vert=1e-3','smag, hori=5e-4,vert=2e-4','Fontsize',fontsize+3)

% legend('smag, hori=1e-4,vert=2e-4','no smag, hori=1e-4,vert=2e-4','no smag, hori=vert=2e-4',...
%     'no smag, hori=2e-4,vert=1e-4','no smag, hori=4e-4,vert=2e-4','no smag, hori=6e-4,vert=3e-4')

