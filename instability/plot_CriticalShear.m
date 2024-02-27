
clear

topo_parm = 1:0.01:20;
Ptide_parm = [1 2]*43200;

Ptide = Ptide_parm(1);
omega = 2*pi/Ptide;

NT = length(topo_parm);
for nt = 1:NT
    topo = topo_parm(nt);    
    critical_convec_shear(nt) = cotd(topo)*omega;
end

figure(10)
plot(topo_parm,critical_convec_shear)
grid on;
grid minor;



%%



clear;
load('CriticalShear.mat')
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


figure(2)
pcolor(log10(Nsquare_parm),topo_parm,log10(criticalShear_M2)');
shading interp;colorbar;
colormap(WhiteBlueGreenYellowRed(0))

figure(3)
clf;
% plot(log10(Bu(10,:)),log10(criticalShear_M2(10,:)))
% hold on;
% plot(log10(Bu(11,:)),log10(criticalShear_M2(11,:)))
% plot(log10(Bu(12,:)),log10(criticalShear_M2(12,:)))
% plot(log10(Bu(13,:)),log10(criticalShear_M2(13,:)))
% plot(log10(Bu(14,:)),log10(criticalShear_M2(14,:)))

% loglog(Bu(50,:),criticalShear_M2(50,:))
%%% (Stratification,slope)
loglog(Bu',criticalShear_M2')

figure(4)
clf;
pcolor(log10(Nsquare_parm),topo_parm,log10(Bu)');
shading interp;colorbar;
colormap(WhiteBlueGreenYellowRed(0))

