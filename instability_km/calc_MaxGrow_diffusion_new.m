
% clear;close all

% expdir = 'exps_new/topo4_kv2e-4/';
% shear_all = (0:0.1:1.8)/1e3;

expdir = 'exps_new/flat_kv2e-4/';
shear_all = (0:0.1:2.0)/1e3;

for ns=1:length(shear_all)
    shear = shear_all(ns);
    load([expdir 'shear' num2str(shear*1e3,3) '_output.mat']);
    [max_grow(ns),I] = max(grow_s,[],'all','omitnan');
    [kk_idx,mm_idx] = ind2sub(size(grow_s),I);
    max_m0(ns) = m0_all(mm_idx);
    max_kx(ns) = kx_all(kk_idx);
    max_lz(ns) = lam_z_all(mm_idx);
    max_lx(ns) = lam_x_all(kk_idx);
    % figure(2)
    % pcolor(kx_all,m0_all,grow_s');
    % shading flat;colorbar;colormap(redblue);
    % clim([0 0.25]);xlim([0 0.06]);ylim([0 0.5])
end

figure(1)
plot(shear_all,max_grow);

% save('../figures/fig4/topo4_km_kv2e-4.mat')
save('../figures/fig4/flat_km_kv2e-4.mat')


%%
figure(2)
subplot(2,2,1)
plot(shear_all,max_lz);xlabel('Shear (1/s)');
title('\lambda_z of the most unstable mode (m)')
set(gca,'Fontsize',16);grid on;grid minor
subplot(2,2,2)
plot(shear_all,max_lx);xlabel('Shear (1/s)')
title('\lambda_x of the most unstable mode (m)')
set(gca,'Fontsize',16);grid on;grid minor
subplot(2,2,3)
plot(shear_all,max_m0./max_kx);xlabel('Shear (1/s)')
title('m_0/k_0 of the most unstable mode')
set(gca,'Fontsize',16);grid on;grid minor


