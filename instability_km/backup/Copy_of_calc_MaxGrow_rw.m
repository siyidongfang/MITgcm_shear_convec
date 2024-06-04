
clear;close all

constants;

expdir = 'parallel_flat_rw/';
shear_all = (0:0.1:2.0)/1e3;

fix_Lx = 3000;
fix_k0 = 2*pi/fix_Lx;
calc_m0_all = fix_k0./rw_all;
calc_lam_x_all = 2*pi./calc_m0_all;
rw_max = 250/fix_Lx;

grow_sr = NaN*zeros(length(shear_all),length(rw_all));

for s = 1:length(shear_all)
    shear = shear_all(s);

    % rw_max_s(s) = 1/(shear/omega+1/rw_max);
    % rw_idx_s = find(rw_all>rw_max_s(s));

    for i=1:Nrw
        fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(i) '.mat'];
        load(fname,'grow');
        grow_sr(s,i)=grow(i);
    end
    % grow_sr(s,rw_idx_s)=NaN;
end


rw_idx = 1:Nrw;
% rw_idx = find(rw_all<=rw_max);
% rw_idx = 75:Nrw;
% rw_idx = 100:Nrw;
figure(1)
pcolor(shear_all,(lam_x_real(rw_idx)),grow_sr(:,rw_idx)');shading flat;colormap(WhiteBlueGreenYellowRed(0))
colorbar;clim([0 0.3])

[max_grow,I]= max(grow_sr(:,rw_idx),[],2);
figure(2)
plot(shear_all,max_grow)

rw_crop = rw_all(rw_idx);
rw_max = rw_crop(I);
% save('grow_rw.mat','max_grow','rw_max')
