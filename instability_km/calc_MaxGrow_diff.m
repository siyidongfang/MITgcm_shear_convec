
clear;

% expdir = 'exps_topo4_diff/';
% shear_all = (0:0.2:1.8)/1e3;

expdir = 'exps_flat_diff/';
shear_all = (0:0.1:1.3)/1e3;

m0_all = [0:0.1:10];
kx_all = [-0.5:0.0025*4:0.5];

for s = 1:length(shear_all)
    shear = shear_all(s)
    for m=1:length(m0_all)
        m0 = m0_all(m);
        load([[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.mat'],'grow')
        grow_smk(s,m,:) = grow;
    end
    figure(1)
    pcolor(kx_all,m0_all,squeeze(grow_smk(s,:,:)));shading flat;colorbar;colormap(redblue);clim([-0.3 0.3])
    max_grow(s) = max(grow_smk(s,:,:),[],'all');
end

figure(2)
plot(shear_all,max_grow);


