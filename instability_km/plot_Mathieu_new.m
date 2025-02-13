%%% Plot growth rate as a function of omega0^2 and epsilon -- as the
%%% classical Arnold instability tongues

clear;addpath ../analysis/colormaps/

expdir = 'exps_new/topo0_nu0_Mathieu'
load([expdir '_output.mat'])

grow_all = [];
omega0_all = [];
epsilon_all = [];
Ri_all = [];
omega0 = N*kx_all/m0_rw;

dt_ri = 10;
Nt_ri = 100*Ptide/dt_ri;
tt_ri = dt_ri:dt_ri:Nt_ri*dt_ri;

Ri_shear = NaN.*zeros(1,Ns); %%% minimum Ri as a function of shear

for ns=1:Ns
    shear = shear_all(ns)

    epsilon = 2*omega0.^3/omega*shear/N;
    epsilon_all=[epsilon_all epsilon];

    grow = grow_rw(ns,:);
    grow_all = [grow_all grow];

    omega0_all = [omega0_all omega0];


    %%% Compute Richardson number
    Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N^2*cosd(topo) - N^2*sind(topo)/omega*shear*sin(omega*tt_ri));
    % figure(50);plot(tt_ri,Ri_inverse)
    Ri_shear(ns) = 1/max(Ri_inverse);  
    if(min(Ri_inverse)<0)
        isConvec = 1
    end

    Ri_all = [Ri_all Ri_shear(ns).*ones(1,Nk)];

end


% for ns = 1:100
%     shear = shear_all(ns);
%     load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01.mat'],'grow','rw_all','m0','N','omega')
%     k0_all = rw_all*m0;
%     omega0 = N*k0_all/m0;
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
%     % ylim([-0.1 3.2]*1e-2)
%     % saveas(fig,[expdir '/figs/topo4_shearN' num2str(ns) '.jpeg']);
% 
% end


[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);
Ri_sort = Ri_all(I);

save('../figures/fig5/fig5_topo0_Ri.mat')



figure(1)
clf;set(gcf,'Color','w');
scatter(omega0_sort/omega,epsilon_sort/omega^2,20,grow_sort,'filled')
colorbar;
clim([-0.3 0.3]/100)
colormap(cmocean('balance'))
ylabel('$\sqrt\epsilon/\omega = \sqrt {2\Lambda\omega_0^3/(N\omega^3)}$','Interpreter','latex')
title('Growth rate (1/hour)')
xlabel('$\sqrt\delta/\omega = \omega_0 /\omega = N k_x/m_z/\omega$','Interpreter','latex')
set(gca,'fontsize',20)
xlim([0 2.5])










