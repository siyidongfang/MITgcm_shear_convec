
clear;
% close all;
% % load('MAVS2_velocity.mat','uu_tilde','time','depth','topo')
% load('MAVS2_velocity.mat')
% zidx = 4:18; %%% Full depth, MAVS 2
% % zidx = 12:18; %%% Bottom 100m, MAVS 2

load('MAVS1_velocity.mat')
zidx = 5:18;%%% For MAVS1, full depth
% zidx = 12:18;%%% For MAVS1, bottom 96m

u1 = uu_tilde*cosd(topo);
u2 = ww_tilde*sind(topo);

fontsize = 20;

ratio = mean(abs(u2),'all','omitnan')/mean(abs(u1),'all','omitnan')
%--- the magnitude of u2 is only 1.08% of u1
% figure(1)
% pcolor(uu(:,300:500));shading flat;colorbar;clim([-0.3 0.3]);colormap(redblue)
% 
% figure(2)
% pcolor(u1(:,300:500));shading flat;colorbar;clim([-0.3 0.3]);colormap(redblue)
% figure(3)
% pcolor(u2(:,300:500));shading flat;colorbar;clim([-0.3 0.3]/80);colormap(redblue)


%%

%%% Calculate the linear fit of uu_tilde for each timestep
dt = diff(time);


for i=1:length(time)
    % i

    % u_fit = uu_tilde(zidx,i);
    u_fit = uu(zidx,i);
    depth_fit = depth(zidx);
    
    % u_fit = uu_tilde(:,i);
    % u_fit = u_fit(~isnan(u_fit));
    % depth_fit = depth(~isnan(uu_tilde(:,i)));

    [p,S] = polyfit(depth_fit,u_fit,1); 
    p1(i) = p(1);p2(i) = p(2);
    
    % figure(1)
    % clf;set(gcf,'color','w')
    % plot(u_fit,depth_fit,'LineWidth',2)
    % hold on;
    % plot(depth_fit*p1(i)+p2(i),depth_fit,'k--','LineWidth',2)
    % hold off;
    % grid on;grid minor;
    % ylabel('Depth (m)');xlabel('u (m/s)');set(gca,'FontSize',fontsize)
    % axis ij;
    % title('Along-canyon velocity')
    % ylim([depth(3) depth(19)]);xlim([-0.3 0.4])

end

shear_linear = -p1;

u_reconstruct = -depth'.*shear_linear+p2;

save('MAVS1_LinearShear.mat','depth','uu','uu_tilde','topo','u_reconstruct','shear_linear','p1','p2','time')
% save('MAVS2_LinearShear_100m.mat','depth','uu_tilde','topo','u_reconstruct','shear_linear','p1','p2','time')


%%

figure(1)
clf;set(gcf,'Color','w');
plot(time/86400,shear_linear);
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Time (days)')
ylabel('shear (1/s)')
title('Linear-fit shear at MAVS2')
ylim([-2.2 1.5]/1e3)
xlim([0 90])


% t_idx_plot = 1:length(time);

t_idx_plot = 1:900;
figure(2);
clf;set(gcf,'color','w','Position',[249 367 931 599]);
subplot(2,1,1)
pcolor(time(t_idx_plot)/86400,depth(zidx),u_reconstruct(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Along-canyon velocity: reconstructed using the linear-fit shear (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;

subplot(2,1,2)
pcolor(time(t_idx_plot)/86400,depth(zidx),uu_tilde(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Along-canyon velocity: observed (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;





