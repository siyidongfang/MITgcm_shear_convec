clear;close all;

load_colors;
fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [100 153 1000 900])
tiledlay = tiledlayout(3,2);

load('grow_M2_calc.mat','isConvec','grow_eig','topo_all','N_all','Ri','Nn','Ntopo','shear_all','Ns')
aaa=isnan(isConvec);
aaa=double(aaa);
aaa(aaa==0)=NaN;
grow_eig = grow_eig.*aaa;

colorn1 = WhiteBlueGreenYellowRed(0);
colorn1 = colorn1(25:25:end,:);
colorn = cmocean('balance');
colorn = colorn(25:30:end,:);
% colorn = colorn1;

nexttile(1);
hold on;
% N=1e-3;
% shear_all = [0:3*N/(Ns-1):3*N];
for n=1:Ntopo
plot(1./squeeze(Ri(n,3,:)),squeeze(grow_eig(n,3,:)),'LineWidth',2,'Color',colorn1(n,:));
% plot(squeeze(shear_all),squeeze(grow_eig(n,3,:)),'LineWidth',2)
end
ylabel('Growth rate (hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('M2 tide, varying topographic slope $\theta$','interpreter','latex','Fontsize',fontsize+5);
box on;grid on;grid minor;
xlim([0 9])
ylim([0 0.805])
colormap(colorn1);
clim([-1 19])
h1=colorbar('XTickLabel',{'0','2','4','6','8','10','12','14','16','18'},'XTick', topo_all);
set(get(h1,'Title'),'String','\ \ \ \ \ \ $(^\circ)$','Fontsize',fontsize,'interpreter','latex');



nexttile(3)
hold on;
for n=1:Nn
plot(1./squeeze(Ri(1,n,:)),squeeze(grow_eig(1,n,:)),'LineWidth',2,'Color',colorn(n,:));
end
ylabel('Growth rate (hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('M2 tide, flat bottom $\theta=0^\circ$, varying $\tilde N$','interpreter','latex','Fontsize',fontsize+5);
box on;grid on;grid minor;
xlim([0 9])
ylim([0 0.805])
colormap(colorn);
clim([0.5 8.5])
h3=colorbar('XTickLabel',{'1e-4','5e-4','1e-3','3e-3','5e-3','7e-3','9e-3','1.1e-2'},'XTick',[1:1:8]);
set(get(h3,'Title'),'String','\ \ \ \ \ \ \ \ $(s^{-1})$','Fontsize',fontsize,'interpreter','latex');




nexttile(5)
hold on;
for n=1:Nn
plot(1./squeeze(Ri(3,n,:)),squeeze(grow_eig(3,n,:)),'LineWidth',2,'Color',colorn(n,:));
end
ylabel('Growth rate (hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('M2 tide, sloping bottom $\theta=4^\circ$, varying $\tilde N$','interpreter','latex','Fontsize',fontsize+5);
box on;grid on;grid minor;
xlim([0 9])
% ylim([0 1.5])
ylim([0 0.805])
colormap(colorn);
clim([0.5 8.5])
h3=colorbar('XTickLabel',{'1e-4','5e-4','1e-3','3e-3','5e-3','7e-3','9e-3','1.1e-2'},'XTick',[1:1:8]);
set(get(h3,'Title'),'String','\ \ \ \ \ \ \ \ \ \ \ \ \ \ $(s^{-1})$','Fontsize',fontsize,'interpreter','latex');




clear isConvec grow_eig topo_all N_all Ri Nn Ntopo
load('grow_K1_calc.mat','isConvec','grow_eig','topo_all','N_all','Ri','Nn','Ntopo')
aaa=isnan(isConvec);
aaa=double(aaa);
aaa(aaa==0)=NaN;
grow_eig = grow_eig.*aaa;


nexttile(2)
hold on;
for n=1:Ntopo
plot(1./squeeze(Ri(n,3,:)),squeeze(grow_eig(n,3,:)),'LineWidth',2,'Color',colorn1(n,:));
end
ylabel('Growth rate (hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('K1 tide, varying topographic slope $\theta$','interpreter','latex','Fontsize',fontsize+5);
box on;grid on;grid minor;
xlim([0 9])
ylim([0 0.805])
colormap(colorn1);
clim([-1 19])
h2=colorbar('XTickLabel',{'0','2','4','6','8','10','12','14','16','18'},'XTick',topo_all);
set(get(h2,'Title'),'String','\ \ \ \ \ \ $(^\circ)$','Fontsize',fontsize,'interpreter','latex');



nexttile(4)
hold on;
for n=1:Nn
plot(1./squeeze(Ri(1,n,:)),squeeze(grow_eig(1,n,:)),'LineWidth',2,'Color',colorn(n,:));
end
ylabel('Growth rate (hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('K1 tide, flat bottom $\theta=0^\circ$, varying $\tilde N$','interpreter','latex','Fontsize',fontsize+5);
box on;grid on;grid minor;
xlim([0 9])
ylim([0 0.805])
colormap(colorn);
clim([0.5 8.5])
h3=colorbar('XTickLabel',{'1e-4','5e-4','1e-3','3e-3','5e-3','7e-3','9e-3','1.1e-2'},'XTick',[1:1:8]);
set(get(h3,'Title'),'String','\ \ \ \ \ \ \ \ $(s^{-1})$','Fontsize',fontsize,'interpreter','latex');





nexttile(6)
hold on;
for n=1:Nn
plot(1./squeeze(Ri(3,n,:)),squeeze(grow_eig(3,n,:)),'LineWidth',2,'Color',colorn(n,:));
end
ylabel('Growth rate (hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
title('K1 tide, sloping bottom $\theta=4^\circ$, varying $\tilde N$','interpreter','latex','Fontsize',fontsize+5);
box on;grid on;grid minor;
xlim([0 9])
ylim([0 0.805])
colormap(colorn);
clim([0.5 8.5])
h3=colorbar('XTickLabel',{'1e-4','5e-4','1e-3','3e-3','5e-3','7e-3','9e-3','1.1e-2'},'XTick',[1:1:8]);
set(get(h3,'Title'),'String','\ \ \ \ \ \ \ \ \ \ \ \ \ \ $(s^{-1})$','Fontsize',fontsize,'interpreter','latex');



% tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal','Direction','TopDown')


% print('-dpng','-r300','fig_supp/figS_sens_eig_matlab2.png');
