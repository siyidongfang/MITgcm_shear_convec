
clear;close all

constants;

grow_sr = NaN*zeros(length(shear_all),length(rw_all));

for s = 1:length(shear_all)
    shear = shear_all(s);

    for i=1:Nrw
        fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(i) '.mat'];
        load(fname,'grow');
        grow_sr(s,i)=grow(i);
    end
end

% fix_Lx = 3000;
% fix_k0 = 2*pi/fix_Lx;
% calc_m0_all = fix_k0./rw_all;
% calc_lam_x_all = 2*pi./calc_m0_all;
% rw_max = 250/fix_Lx;

% rw_idx = find(rw_all<=rw_max);

% nstart = 84+61;
nstart = 84+61+23;
% nstart=1
rw_idx = nstart:Nrw;
lam_x_real(nstart)
figure(1)
set(gcf,'Color','w');pcolor(shear_all,(lam_x_real(rw_idx))/1000,grow_sr(:,rw_idx)');set(gca,'Fontsize',20);shading flat;colormap(WhiteBlueGreenYellowRed(0))
colorbar;clim([0 0.35]);xlabel('Shear (1/s)');ylabel('Horizontal wavelength (km)');title('Growth rate (1/hour)','Fontsize',25)
% ylim([0 30000])

[max_grow,I]= max(grow_sr(:,rw_idx),[],2);
figure(2)
set(gcf,'Color','w');plot(shear_all,max_grow,'LineWidth',2);set(gca,'Fontsize',20);grid on;xlabel('Shear (1/s)'); ylabel('Growth rate (1/hour)')

rw_crop = rw_all(rw_idx);
rw_max = rw_crop(I);

save('grow_rw_new_3km.mat')
% save('grow_rw_new_5.6km.mat')
% save('grow_rw_new.mat')

