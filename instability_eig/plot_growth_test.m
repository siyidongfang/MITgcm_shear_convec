clear;
expdir = 'test1/';

Ptide = 43200;
N = sqrt(1)*1e-3;
topo=0;

h_shear = 250;
lam_z_all = h_shear;
lam_x_all = [0:20:9000];
Nlam = length(lam_x_all);
shear_all = [0.5e-3:0.025e-3:1.6e-3]; 
Ns = length(shear_all);

r = lam_x_all/h_shear;

shear_floquet = shear_all;
grow_all = zeros(Ns,Nlam);
for ns=1:Ns
    shear = shear_all(ns);
    load([expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3) '_shear' num2str(shear*1e3) '.mat'],'grow','grow_pos','grow_neg');

    grow(grow<=0)=NaN;
    grow = log(grow)/43200*3600;

    % grow_pos(grow_pos<=0)=NaN;
    % grow_pos = log(grow_pos)/43200*3600;
    % 
    % grow_neg(grow_neg<=0)=NaN;
    % grow_neg = log(grow_neg)/43200*3600;

    % grow_pos_all(ns,:) = grow_pos;
    % grow_neg_all(ns,:) = grow_neg;
    grow_all(ns,:) = grow;

    max_grow_floquet(ns) =max(grow,[],'all','omitnan');
end

figure(1)
clf;
set(gcf,'Color','w')

% subplot(1,3,1)
% pcolor(shear_all,r,grow_pos_all');shading flat;colorbar;
% colormap(redblue);
% ylim([0 30])
% clim([-5 5]/10)
% 
% subplot(1,3,2)
% pcolor(shear_all,r,grow_neg_all');shading flat;colorbar;
% colormap(redblue);
% ylim([0 30])
% clim([-5 5]/10)


subplot(1,3,3)
pcolor(shear_all,r,grow_all');shading flat;colorbar;
colormap(redblue);
ylim([0 30])
clim([-5 5]/10)



% figure(2)
% clf;
% set(gcf,'Color','w')
% plot(shear_floquet,max_grow_floquet,'LineWidth',2)
% set(gca,'fontsize',17)
% xlabel('shear (1/s)')
% ylabel('(1/hour)')
% title('The Floquet exponents as a function of shear')
% ylim([0 0.34])
% grid on;grid minor;




