
clear;close all

constants;
% load_colors;

expdir = 'parallel_flat_rw/';
shear_all = (0:0.1:2.0)/1e3;

grow_sr = NaN*zeros(length(shear_all),length(rw_all));


for s = 1:length(shear_all)
    shear = shear_all(s)
    for i=1:Nrw
        fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_rw' num2str(i) '.mat'];
        load(fname,'grow');
        grow_sr(s,i)=grow(i);
    end
end

figure(1)
pcolor(shear_all,(lam_x_real),grow_sr');shading flat;colormap(WhiteBlueGreenYellowRed(0))

max_grow = max(grow_sr,[],2);
figure(2)
plot(shear_all,max_grow)
