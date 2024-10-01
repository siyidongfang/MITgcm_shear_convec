load('GrowthRate_new_flat_tanh.mat')

figure(1)
subplot(1,2,1)
pcolor(GrowthRate_Floquet)
shading flat;colorbar;colormap(redblue);

subplot(1,2,2)
plot(shear_Floquet,growth_Floquet)
