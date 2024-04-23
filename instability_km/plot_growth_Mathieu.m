%%% Plot growth rate as a function of omega0^2 and epsilon -- as the
%%% classical Arnold instability tongues

clear;addpath ../analysis/colormaps/

%%% calculate omega0 and epsilon of each case
% expdir = 'exps_flatbottom';
expdir = 'exps_flat';
shear_all = [0:0.005:0.35 0.5:0.005:0.6 0.61:0.005:0.985 1:0.005:2]*1e-3; % 1:0.005:1.5
% shear_all = [0:0.005:0.32]*1e-3;
grow_all = [];
omega0_all = [];
epsilon_all = [];

for ns = 1:92
    shear = shear_all(ns);
    load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01.mat'],'grow','rw_all','m0','N','omega')
    k0_all = rw_all*m0;
    omega0 = N*k0_all*m0;
    epsilon = 2*omega0.^2.*omega0/omega*shear/N;
    epsilon = (epsilon/omega^2).^(1/2);
    grow_all = [grow_all grow];
    omega0_all = [omega0_all omega0];
    epsilon_all = [epsilon_all epsilon];

    fig=figure(2);
    set(fig,'visible','on');
    clf;set(gcf,'Color','w');
    plot((rw_all)*N/omega,grow,'LineWidth',2)
    xlabel('Natural frequency/Forcing frequency (N k_x/m_z/\omega)')
    title('Growth rate (1/hour)')
    grid on;grid minor;
    ylabel('(1/hour)')
    set(gca,'fontsize',20)
    set(gcf, 'InvertHardcopy', 'off')
    ylim([-0.1 10]*1e-2)
    saveas(fig,[expdir '/figs/flat_shearN' num2str(ns) '.jpeg']);

end



% for ns =169:217
%     shear = shear_all(ns);
%     load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01.mat'],'grow','rw_all','m0','N','omega')
%     k0_all = rw_all*m0;
%     omega0 = N*k0_all*m0;
%     epsilon = 2*omega0.^2.*omega0/omega*shear/N;
%     epsilon = (epsilon/omega^2).^(1/2);
%     grow_all = [grow_all grow];
%     omega0_all = [omega0_all omega0];
%     epsilon_all = [epsilon_all epsilon];
% 
%     % fig=figure(2);
%     % set(fig,'visible','on');
%     % clf;set(gcf,'Color','w');
%     % plot((rw_all)*N/omega,grow,'LineWidth',2)
%     % xlabel('Natural frequency/Forcing frequency (N k_x/m_z/\omega)')
%     % title('Growth rate (1/hour)')
%     % grid on;grid minor;
%     % ylabel('(1/hour)')
%     % set(gca,'fontsize',20)
%     % set(gcf, 'InvertHardcopy', 'off')
%     % ylim([-0.1 10]*1e-2)
%     % saveas(fig,[expdir '/figs/flat_shearN' num2str(ns+22) '.jpeg']);
% 
% end
% 
% for ns = 218:251
%     shear = shear_all(ns);
%     load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01_newRW.mat'],'grow','rw_all','m0','N','omega')
%     k0_all = rw_all*m0;
%     omega0 = N*k0_all*m0;
%     epsilon = 2*omega0.^2.*omega0/omega*shear/N;
%     epsilon = (epsilon/omega^2).^(1/2);
%     grow_all = [grow_all grow];
%     omega0_all = [omega0_all omega0];
%     epsilon_all = [epsilon_all epsilon];
% 
%     % fig=figure(2);
%     % set(fig,'visible','on');
%     % clf;set(gcf,'Color','w');
%     % plot((rw_all)*N/omega,grow,'LineWidth',2)
%     % xlabel('Natural frequency/Forcing frequency (N k_x/m_z/\omega)')
%     % title('Growth rate (1/hour)')
%     % grid on;grid minor;
%     % ylabel('(1/hour)')
%     % set(gca,'fontsize',20)
%     % set(gcf, 'InvertHardcopy', 'off')
%     % ylim([-0.1 10]*1e-2)
%     % saveas(fig,[expdir '/figs/flat_shearN' num2str(ns+22) '.jpeg']);
% 
% end
% 
% for ns = 252:length(shear_all)
%     shear = shear_all(ns);
%     load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01_newRW2.mat'],'grow','rw_all','m0','N','omega')
%     k0_all = rw_all*m0;
%     omega0 = N*k0_all*m0;
%     epsilon = 2*omega0.^2.*omega0/omega*shear/N;
%     epsilon = (epsilon/omega^2).^(1/2);
%     grow_all = [grow_all grow];
%     omega0_all = [omega0_all omega0];
%     epsilon_all = [epsilon_all epsilon];
% 
%     fig=figure(2);
%     set(fig,'visible','on');
%     clf;set(gcf,'Color','w');
%     plot((rw_all)*N/omega,grow,'LineWidth',2)
%     xlabel('Natural frequency/Forcing frequency (N k_x/m_z/\omega)')
%     title('Growth rate (1/hour)')
%     grid on;grid minor;
%     ylabel('(1/hour)')
%     set(gca,'fontsize',20)
%     set(gcf, 'InvertHardcopy', 'off')
%     ylim([-0.1 40]*1e-2)
%     saveas(fig,[expdir '/figs/flat_shearN' num2str(ns+22) '.jpeg']);
% 
% end

[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);

figure(1)
clf;set(gcf,'Color','w');
scatter(omega0_sort/omega,epsilon_sort/omega^2,20,grow_sort,'filled')
colorbar;
% xlim([0 1.5])
% ylim([0 0.5])
% clim([-0.3 0.3]/10)
colormap(cmocean('balance'))
% colormap(redblue)
ylabel('$\sqrt\epsilon/\omega = \sqrt {2\Lambda\omega_0^3/(N\omega^3)}$','Interpreter','latex')
title('Growth rate (1/hour)')
xlabel('$\sqrt\delta/\omega = \omega_0 /\omega = N k_x/m_z/\omega$','Interpreter','latex')
set(gca,'fontsize',20)











