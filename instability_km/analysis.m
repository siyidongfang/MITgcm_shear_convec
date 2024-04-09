
% clear;
% close all;
% shear_all = [0:1e-4:1.8e-3]; % Ri=1, shear = 0.97e-3;
% rw_all= 10.^([-2:0.1:-1 -0.95:0.01:-0.5 0.6:0.1:1]);
shear_all = [0:0.1:1.8]*1e-3;
% growrate = zeros(length(shear_all),length(rw_all));
for ns = 1:length(shear_all)
    shear = shear_all(ns)
    load(['output/topo4_Nsq1e-6/growth_shear' num2str(shear*1e3,3) '.mat'])
    length_grow = length(grow);
    growrate(ns,1:length_grow)=grow;
    [max_growth(ns) rw_idx] = max(grow);
    rw_mg(ns) = rw_all(rw_idx);
end

shear_all = [0:0.1:1.8]*1e-3;

%%% Find out the wavenumber ratio rw=kx/mz corresponding to the maximum
%%% growth rate, for each shear value

% growth_round = growrate;
% growth_round(growrate>0.2) = round(growrate(growrate>0.2),1);


% figure(22);
% clf;
% set(gcf,'Color','w');
% % pcolor(shear_all,atand(1./rw_all),growrate')
% pcolor(shear_all,log10(rw_all),growrate')
% shading flat;colormap(WhiteBlueGreenYellowRed(0))
% hold on;
% scatter(shear_all,log10(rw_mg),50,'filled','black')
% grid on;grid minor;
% title('Growth rate (1/hour)')
% xlabel('Shear (1/s)')
% ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
% set(gca,'fontsize',20)
% colorbar;
% clim([0 0.25])

figure(23)
set(gcf,'Color','w')
l3 = plot(shear_all,max_growth,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
hold on;
plot(shear_all,max_growth*0,'k--');
xlabel('Shear (1/s)')
title('Maximum growth rate (1/hour)')
ylabel('(1/hour)')
legend([l1 l2 l3],'$\nu=\kappa=0$',...
    '$\nu=2\times 10^{-6}m^2/s,\,\kappa=1\times 10^{-6}m^2/s$',...
    '$\nu=1\times 10^{-5}m^2/s,\,\kappa=1\times 10^{-6}m^2/s$',...
    'Interpreter','latex','Fontsize',30)

