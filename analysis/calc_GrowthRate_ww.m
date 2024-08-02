%%% 
%%% calc_GrowthRate_RMSE.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations


clear;

for  ne = 1:4
load_all

No = 40*12;
tidx = 1:No;

fit_span = 10*12+1:No;
% fit_span = 3*12+1:9*12;

Nt = length(tidx);
Hshear = 250;
dz = delR(end);
zidx = 300:500;
Nshear = length(zidx);
div_ww = zeros(Nt,Nshear);

parfor o = tidx
    nIter = dumpIters(o);
    time_h(o) = nIter.*deltaT./3600;
    time_tidal(o) = time_h(o)/12;

    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    ww_shear = ww(:,zidx);
    ww2(o,:) = sum(dx*ww_shear.^2);
end

ww2_zint = sum(dz*ww2,2);
ke = ww2_zint/2;

xxplot = time_h;
yyplot = log(ke)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow(ne) = pp(1)
[y_fit,delta_fit] = polyval(pp,xxplot,S);

filename = [expdir expname '/RMSE_ww.mat'];
save(filename,'time_h','xxplot','yyplot','fit_span','pp','y_fit','ww2_zint','grow')

h_figure = figure(1);
clf;
set(gcf,'Color','w','Position',[211 289 852 394])
plot(time_h/12,log(ke)/2,'LineWidth',2);
hold on;
plot(xxplot(fit_span)/12, y_fit(fit_span),'k-.','LineWidth',2);
set(gca,'Fontsize',fontsize+5)
xlabel('Time (tidal cycles)','Interpreter','latex')
title('Normalized vertical turbulent energy in the shear layer','Interpreter','latex')
ylabel('$\log(\mathrm{turbulent\ energy})$','Interpreter','latex')
grid on;grid minor;
hold on;

ylim([-40 5])
print('-dpng','-r150',[expdir expname '_ww2.png']);


end
