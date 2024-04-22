%%% Plot growth rate as a function of omega0^2 and epsilon -- as the
%%% classical Arnold instability tongues

clear;addpath ../analysis/colormaps/

%%% calculate omega0 and epsilon of each case
% expdir = 'exps_flatbottom';
expdir = 'exps_test';
shear_all = [0:0.01:1.8]*1e-3;

grow_all = [];
omega0_all = [];
epsilon_all = [];

% for ns = 1:length(shear_all)
for ns = 1:50
    shear = shear_all(ns);
    load([expdir '/growth_shear' num2str(shear*1e3,3) '.mat'],'grow','rw_all','m0','N','omega','topo')
    % load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01_test.mat'],'grow','rw_all','m0','N','omega')
    k0_all = rw_all*m0;
    omega0 = N*k0_all*m0;
    epsilon = 2*omega0.^2.*omega0/omega*shear/N;
    grow_all = [grow_all grow];
    omega0_all = [omega0_all omega0];
    epsilon_all = [epsilon_all epsilon];

    % figure(1)
    % % plot(epsilon/omega^2,grow)
    % plot(omega0.^2/omega^2,grow)
    % hold on;
    % xlim([0 2])
end


[omega0_sort,I] = sort(omega0_all);
epsilon_sort = epsilon_all(I);
grow_sort = grow_all(I);

figure(1)
scatter(omega0_sort.^2/omega^2,epsilon_sort/omega^2,[],grow_sort)
colorbar;
xlim([0 2])
clim([-0.04 0.04])
colormap(cmocean('balance'))








