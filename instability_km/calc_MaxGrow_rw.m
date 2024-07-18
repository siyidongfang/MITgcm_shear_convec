
clear;close all
% constants;
% expdir = 'exps_new/topo0_nu0/'
expdir = 'exps_new/topo0_nu2e-4_lores_'
load([expdir 'output.mat']);

% load([expdir 'grow_rw.mat'])
% grow_sr = NaN*zeros(length(shear_all),length(rw_all));

% for s = 1:length(shear_all)
%     s
%     shear = shear_all(s);
% 
%     for i=1:Nrw
%         fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(i) '.mat'];
%         load_func(fname);
%         grow_sr(s,i)=grow(i);
%     end
% end

% save([expdir '_grow.mat']);

%%

% load([expdir 'grow.mat']);



% h_reset = 1500;
h_reset = 250;
L_reset = h_reset./rw_all;

rw_idx = 1:Nrw;

crop_limit = 3000;
rw_idx=find(L_reset<=crop_limit);



figure(1)
clf;
set(gcf,'Color','w');
pcolor(shear_all,h_reset/1000./(rw_all(rw_idx)),grow_rw(:,rw_idx)');
set(gca,'Fontsize',20);
shading flat;colormap(WhiteBlueGreenYellowRed(0))
% set(gcf,'Color','w');pcolor(shear_all,(lam_x_real(rw_idx))/1000,grow_sr(:,rw_idx)');set(gca,'Fontsize',20);shading flat;colormap(WhiteBlueGreenYellowRed(0))
colorbar;clim([0 0.35]);xlabel('Shear (1/s)');ylabel('Horizontal wavelength (km)');title('Growth rate (1/hour)','Fontsize',25)
% ylim([0 13])
% ylim([0 35])
ylim([0 6])


figure(3)
set(gcf,'Color','w');pcolor(shear_all,kx_all(rw_idx),grow_rw(:,rw_idx)');set(gca,'Fontsize',20);shading flat;colormap(WhiteBlueGreenYellowRed(0))
colorbar;clim([0 0.35]);xlabel('Shear (1/s)');ylabel('Horizontal wavenumber k_0 (1/m)');title('Growth rate (1/hour)','Fontsize',25)
ylim([0 0.07]);
box on;

[max_grow,I]= max(grow_rw(:,rw_idx),[],2);
figure(2)
set(gcf,'Color','w');
% plot(shear_all,max_grow,'LineWidth',2);
scatter(shear_all,max_grow);
set(gca,'Fontsize',20);grid on;xlabel('Shear (1/s)'); ylabel('Growth rate (1/hour)')

rw_crop = rw_all(rw_idx);
rw_max = rw_crop(I);

% save('grow_rw_new_3km.mat')
% save('../figures/fig3/topo4_kappa0.mat')
% save('../figures/fig3/topo4_nu2e-4.mat')



