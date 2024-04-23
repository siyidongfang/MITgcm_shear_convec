%%% Plot growth rate as a function of omega0^2 and epsilon -- as the
%%% classical Arnold instability tongues

clear;addpath ../analysis/colormaps/

%%% calculate omega0 and epsilon of each case
% expdir = 'exps_flatbottom';
expdir = 'exps_Mathieu';
shear_all = [0:0.05:1]*1e-3;

grow_all = [];
omega0_all = [];
epsilon_all = [];

for ns = 1:length(shear_all)
    shear = shear_all(ns);
    load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01.mat'],'grow','rw_all','m0','N','omega')
    k0_all = rw_all*m0;
    omega0 = N*k0_all*m0;
    epsilon = 2*omega0.^2.*omega0/omega*shear/N;
    grow_all = [grow_all grow];
    omega0_all = [omega0_all omega0];
    epsilon_all = [epsilon_all epsilon];

end


[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);

figure(1)
clf;
scatter(omega0_sort.^2/omega^2,epsilon_sort/omega^2,[],grow_sort)
colorbar;
xlim([0 2])
clim([-0.3 0.3]/10)
colormap(cmocean('balance'))










