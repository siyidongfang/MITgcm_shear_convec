
% %%
% fontsize = 20;
% figure(1)
% set(gcf,'color','w','Position',[133 306 717 582]);
% semilogy(Shear_parm,Ri_min,'LineWidth',2)
% hold on;plot(Shear_parm,0.25*ones(1,length(Shear_parm)),'-.','LineWidth',2)
% grid on;
% set(gca,'Fontsize',fontsize);
% ylabel('Ri_{min}')
% xlabel('Shear (1/s)')
% title('Minimum Ri of the background flow','Fontsize',fontsize+3);

%%
figure(2)
pcolor(log10(Nsquare_parm),topo_parm,log10(criticalShear_M2)');
shading flat;colorbar;
colormap(WhiteBlueGreenYellowRed(0))

figure(3)
plot(log10(Bu(10,:)),log10(criticalShear_M2(10,:)))
hold on;
plot(log10(Bu(11,:)),log10(criticalShear_M2(11,:)))
plot(log10(Bu(12,:)),log10(criticalShear_M2(12,:)))
plot(log10(Bu(13,:)),log10(criticalShear_M2(13,:)))
plot(log10(Bu(14,:)),log10(criticalShear_M2(14,:)))


loglog(Bu,criticalShear_M2)
