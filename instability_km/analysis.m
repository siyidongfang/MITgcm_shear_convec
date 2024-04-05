
clear;close all;
shear_all = [0:1e-4:1.8e-3]; % Ri=1, shear = 0.97e-3;
rw_all= 10.^([-2:0.1:-1 -0.95:0.01:-0.5 0.6:0.1:1]);
growrate = zeros(length(shear_all),length(rw_all));
for ns = 1:length(shear_all)
    shear = shear_all(ns)
    load(['output/topo4_Nsq1e-6/growth_shear' num2str(shear*1e3,3) '.mat'])
    growrate(ns,:)=grow;
end

%%% Find out the wavenumber ratio rw=kx/mz corresponding to the maximum
%%% growth rate, for each shear value

% growth_round = growrate;
% growth_round(growrate>0.2) = round(growrate(growrate>0.2),1);
[max_growth rw_idx] = max(growrate,[],2);
rw_mg = rw_all(rw_idx);


figure(22);
clf;
set(gcf,'Color','w');
% pcolor(shear_all,atand(1./rw_all),growrate')
pcolor(shear_all,log10(rw_all),growrate')
shading flat;colormap(WhiteBlueGreenYellowRed(0))
hold on;
scatter(shear_all,log10(rw_mg),50,'filled','black')
grid on;grid minor;
title('Growth rate (1/hour)')
xlabel('Shear (1/s)')
ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
set(gca,'fontsize',20)
colorbar;
% clim([0 6])

figure(23)
clf;set(gcf,'Color','w')
plot(shear_all,max_growth,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
xlabel('Shear (1/s)')
title('Maximum growth rate (1/hour)')
ylabel('(1/hour)')