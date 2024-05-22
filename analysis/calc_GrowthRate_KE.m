%%% 
%%% calc_GrowthRate_KE.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations

clear;
% close all
useTanhShear = false;
useLinearShear = true;
shear_all = (0:0.1:2)*1e-3;

for  ne = 15
load_all

Shear = shear_all(ne);

No = nDumps;
% No = 31*12;
No = 35*12;
tidx = 1:No;
Nt = length(tidx);
h_shear = 250;
dz = delR(end);
Nshear = round(h_shear/dz);
zidx = 300:500;
Nshear = length(zidx);

calc_vrelax;

for o = tidx
    nIter = dumpIters(o);
    time_h(o) = nIter.*deltaT./3600;
    time_tidal(o) = time_h(o)/12;

    t2 = nIter*deltaT;
    if(o>1)
        t1 = dumpIters(o-1)*deltaT;
    else
        t1 = 0;
    end
    utide = 1/(t2-t1)*vrelax/omega*(sin(omega*t2)-sin(omega*t1));

    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));

    usq(o) = sum((uu(:,zidx)-repmat(utide(zidx),[Nx 1])).^2,'all');
    wsq(o) = sum(ww(:,zidx).^2,'all');
    tsq(o) = sum(tt(:,zidx).^2,'all');

end


figure(10)
plot(time_h/12,log10(usq)/2,'LineWidth',2);
hold on;
plot(time_h/12,log10(wsq)/2,'LineWidth',2);
plot(time_h/12,log10(tsq)/2,'LineWidth',2);

%%
Pr = 1;
% ke = usq/2+wsq/2;
ke = wsq/2;
pe = Pr*tsq/2;
energy =ke+pe;

% fit_span = 31*12+1:36*12;
fit_span = 20*12+1:No;

xxplot = time_h;
yyplot = log(energy)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(ne) = pp(1)
[y_fit,delta_fit] = polyval(pp,xxplot,S);


% % filename = [expdir expname '/grow_KE.mat'];
% % save(filename,'time_h','xxplot','yyplot','fit_span','pp','y_fit',...
% %     'energy','ke','pe','div_uu2_zavg','div_vv2_zavg','div_ww2_zavg','div_tt2_zavg','grow')


h_figure = figure(1);
% set(h_figure,'Visible',false);
clf;
set(gcf,'Color','w','Position',[211 289 852 394])
plot(time_h/12,log(pe)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(ke)/2,'LineWidth',2);
plot(time_h/12,log(energy)/2,'LineWidth',2);
plot(xxplot(fit_span)/12, y_fit(fit_span),'--','LineWidth',2);
set(gca,'Fontsize',fontsize)
xlabel('Time (tidal cycles)')
title('Turbulent energy averaged over the bottom shear layer')
legend('RMSE of PE','RMSE of KE','RMSE of PE+KE','Linear fit','Position',[0.6972 0.1497 0.1890 0.2069])
grid on;grid minor;
hold on;
% print('-dpng','-r150',[expdir expname '_grow_ke.png']);


end




