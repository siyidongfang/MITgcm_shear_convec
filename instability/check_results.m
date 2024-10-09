% clear;
% load('products/GrowthRate_new_flat_linear.mat')
% load('products/GrowthRate_new_flat_tanh.mat')
% load('products/GrowthRate_new_topo4_linear.mat')
% load('products/GrowthRate_new_topo4_tanh.mat')
% load('products/GrowthRate_new_topo4_0Center.mat')

figure(1)
clf
subplot(1,2,1)
pcolor(shear_Floquet,log10(lambda_Floquet),GrowthRate_Floquet);colorbar;
subplot(1,2,2)
plot(shear_Floquet,growth_Floquet)
grid on;grid minor;
hold on;
load('../instability_km/exps_varyingN/N1e-3output.mat','shear_all','grow_rw')
[max_grow_nodiff r_idx]=max(grow_rw,[],2);
plot(shear_all,max_grow_nodiff)

xlim([0 2e-3])
