%%%
%%% calc_N2total.m
%%%
%%% Calculate the total stratification

close all;
% load('/Users/ysi/MITgcm_shear_convec/instability/exps_mechanism/lambda100/H300_topo4_Pt43200_N0.001_S0.0013_lambda100/output.mat')
% load('/Users/ysi/MITgcm_shear_convec/instability/exps_mechanism/lambda100/H300_topo4_Pt43200_N0.001_S0.0013_lambda100_reduced/output.mat')
load('/Users/ysi/MITgcm_shear_convec/instability/exps_mechanism/lambda100/H300_topo4_Pt43200_N0.001_S0.0013_lambda100_reduced_NOwN2/output.mat')

% load('/Users/ysi/MITgcm_shear_convec/instability/exps_mechanism/lambda398/H300_topo4_Pt43200_N0.001_S0.0013_lambda398/output.mat')
% load('/Users/ysi/MITgcm_shear_convec/instability/exps_mechanism/lambda398/H300_topo4_Pt43200_N0.001_S0.0013_lambda398_reduced/output.mat')

% pidx = round(Nt/20*3):round(Nt/20*8);
pidx = 1:Nt;
plot_time = ttd(pidx)/43200;
fontsize = 20;

CLIM1 = [-3 3]/1e6;
CLIM2 = CLIM1/1e5;

dtd = dt/omega;
dzd = dz*delta;
zz_dbdz = 0.5*(zzd(1:end-1)+zzd(2:end));
re_buoyd = re_buoy*N^2*U1*sind(topo)/omega;

dbdz = diff(re_buoyd,1,2)/dzd;
dB0_dz = N^2 * cosd(topo) * ones(1,length(ttd));
dB0_dz = repmat(dB0_dz',[1 length(zz_dbdz)]);
dB_dz  = N^2 * (-sind(topo)*Shear/omega*sin(omega*ttd));
dB_dz = repmat(dB_dz',[1 length(zz_dbdz)]);
N2_total = dB0_dz+dB_dz+dbdz;

figure(1);set(gcf,'Color','w','Position',[131 138 1452 524]);
subplot(2,2,1);pcolor(plot_time,zz_dbdz,dB0_dz(pidx,:)');title('$dB_0(x,z)/dz$','Interpreter','latex');shading flat;clim(CLIM1);colormap(redblue);colorbar;set(gca,'Fontsize',fontsize);ylabel('HAB (m)');
subplot(2,2,2);pcolor(plot_time,zz_dbdz,dB_dz(pidx,:)');title('$dB(z,t)/dz$','Interpreter','latex');shading flat;clim(CLIM1);colormap(redblue);colorbar;set(gca,'Fontsize',fontsize);ylabel('HAB (m)');
subplot(2,2,3);pcolor(plot_time,zz_dbdz,dbdz(pidx,:)');title('$db^\prime/dz$','Interpreter','latex');shading flat;clim(CLIM1);colormap(redblue);colorbar;set(gca,'Fontsize',fontsize);xlabel('Time (tidal cycles)');ylabel('HAB (m)');
subplot(2,2,4);pcolor(plot_time,zz_dbdz,N2_total(pidx,:)');title('$N^2_\mathrm{total} = dB_0/dz+dB/dz+db^\prime/dz$','Interpreter','latex');shading flat;clim(CLIM1);colormap(redblue);colorbar;set(gca,'Fontsize',fontsize);xlabel('Time (tidal cycles)');ylabel('HAB (m)');

figure(2);set(gcf,'Color','w','Position',[131 376 865 286]);pcolor(plot_time,zzd_wgrid,www(pidx,:)');title('Vertical velocity, w (m/s)','Interpreter','latex');shading flat;colormap redblue; colorbar;clim([-1 1]/500/10);set(gca,'Fontsize',fontsize);xlabel('Time (tidal cycles)');ylabel('HAB (m)');

% Utided = cos(omega*ttd)'.*Atide; % figure(3);pcolor(ttd/3600,zzd,Utided');shading flat;colormap redblue; colorbar;

% bt1 = -www.*dB_dz';
% bt2 = -www.*N^2*cosd(topo);
% % bt3 = -Utided.*db_dx;

bq_w1dt = real(bq2)*dt*N^2*U1*sind(topo)/omega;
bq_w2dt = real(bq4)*dt*N^2*U1*sind(topo)/omega;
bq_udt = real(bq1)*dt*N^2*U1*sind(topo)/omega;

bq_3dt = real(bq3)*dt*N^2*U1*sind(topo)/omega;
bq_5dt = real(bq5)*dt*N^2*U1*sind(topo)/omega;

% figure();clf;pcolor(plot_time,zzd,bq_w1dt(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e7)
% figure();clf;pcolor(plot_time,zzd,bq_w2dt(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e7)
% figure();clf;pcolor(plot_time,zzd,bq_w1dt(pidx,:)'+bq_w2dt(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e7)
% figure();clf;pcolor(plot_time,zzd,bq_udt(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e7)
% 
% figure();clf;pcolor(plot_time,zzd,bq_3dt(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e7)
% figure();clf;pcolor(plot_time,zzd,bq_5dt(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e7)


cumsumbqw1 = cumsum(bq_w1dt);
cumsumbqw2 = cumsum(bq_w2dt);
cumsumbu = cumsum(bq_udt);
cumsumb3 = cumsum(bq_3dt);
cumsumb5 = cumsum(bq_5dt);

% figure();clf;pcolor(plot_time,zzd,cumsumbqw1(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e5)
% figure();clf;pcolor(plot_time,zzd,cumsumbqw2(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e5)
% figure();clf;pcolor(plot_time,zzd,cumsumbqw1(pidx,:)'+cumsumbqw2(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e5)
% figure();clf;pcolor(plot_time,zzd,cumsumbu(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e5)
% figure();clf;pcolor(plot_time,zzd,cumsumb3(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e5)
% figure();clf;pcolor(plot_time,zzd,cumsumb5(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/1e5)

db_dz_w1 = diff(cumsumbqw1,1,2)/dzd;
db_dz_w2 = diff(cumsumbqw2,1,2)/dzd;
db_dz_u = diff(cumsumbu,1,2)/dzd;
db_dz_3 = diff(cumsumb3,1,2)/dzd;
db_dz_5 = diff(cumsumb5,1,2)/dzd;

db_dz_reconst = db_dz_w1+db_dz_w2+db_dz_u + db_dz_3+db_dz_5;

figure(3);set(gcf,'Color','w','Position',[131 138 1484 757]);
subplot(3,2,1);pcolor(plot_time,zz_dbdz,db_dz_w1(pidx,:)');title('$\partial_z \int (-w^\prime \partial_z B_0 ) dt = \partial_z \int (-w^\prime \tilde N^2 \cos\theta ) dt$','Interpreter','latex');shading flat;colormap redblue; colorbar;clim(CLIM2);set(gca,'Fontsize',fontsize);ylabel('HAB (m)');
subplot(3,2,2);pcolor(plot_time,zz_dbdz,db_dz_w2(pidx,:)');title('$\partial_z \int (-w^\prime \partial_z B ) dt = \partial_z \int [-w^\prime \tilde N^2 \Lambda/\omega\sin(\omega t) \sin\theta ] dt $','Interpreter','latex');shading flat;colormap redblue; colorbar;clim(CLIM2);set(gca,'Fontsize',fontsize);ylabel('HAB (m)');
subplot(3,2,3);pcolor(plot_time,zz_dbdz,db_dz_u(pidx,:)');title('$\partial_z \int (-U b^\prime_x ) dt $','Interpreter','latex');shading flat;colormap redblue; colorbar;clim(CLIM2);set(gca,'Fontsize',fontsize);ylabel('HAB (m)');
subplot(3,2,4);pcolor(plot_time,zz_dbdz,db_dz_3(pidx,:)');title('$\partial_z \int (-u^\prime \partial_x B_0 ) dt = \partial_z \int (-u^\prime \tilde N^2 \sin\theta ) dt $','Interpreter','latex');shading flat;colormap redblue; colorbar;clim(CLIM2);set(gca,'Fontsize',fontsize);ylabel('HAB (m)');
subplot(3,2,5);pcolor(plot_time,zz_dbdz,db_dz_5(pidx,:)');title('$\partial_z \int (\kappa\nabla^2 b^\prime ) dt $','Interpreter','latex');shading flat;colormap redblue; colorbar;clim(CLIM2);set(gca,'Fontsize',fontsize);xlabel('Time (tidal cycles)');ylabel('HAB (m)');
subplot(3,2,6);pcolor(plot_time,zz_dbdz,db_dz_reconst(pidx,:)');title('$db^\prime/dz$','Interpreter','latex');shading flat;colormap redblue; colorbar;clim(CLIM2);set(gca,'Fontsize',fontsize);xlabel('Time (tidal cycles)');ylabel('HAB (m)');



N2_total_reconst = dB0_dz+dB_dz+db_dz_reconst;

% figure(17);pcolor(plot_time,zz_dbdz,N2_total_reconst(pidx,:)');shading flat;clim([-3 3]/1e6);colormap(redblue);colorbar


% % % % %%
% % % % % db_dz_w1_zavg = -mean(db_dz_w1(db_dz_w1<0),2);
% % % % % db_dz_w2_zavg = -mean(db_dz_w2(db_dz_w2<0),2);
% % % % % db_dz_u_zavg = -mean(db_dz_u(db_dz_u<0),2);
% % % % % db_dz_3_zavg = -mean(db_dz_3(db_dz_3<0),2);
% % % % % db_dz_5_zavg = -mean(db_dz_5(db_dz_5<0),2);
% % % % % db_dz_reconst_zavg = -mean(db_dz_reconst(db_dz_reconst<0),2);
% % % % % dbdz_zavg = -mean(dbdz(dbdz<0),2);
% % % % 
% % % % 
% % % % % db_dz_w1_zavg = mean(db_dz_w1,2);
% % % % % db_dz_w2_zavg = mean(db_dz_w2,2);
% % % % % db_dz_u_zavg = mean(db_dz_u,2);
% % % % % db_dz_3_zavg = mean(db_dz_3,2);
% % % % % db_dz_5_zavg = mean(db_dz_5,2);
% % % % % db_dz_reconst_zavg = mean(db_dz_reconst,2);
% % % % % dbdz_zavg = mean(dbdz,2);
% % % % 
% % % % zzidx = 1;
% % % % 
% % % % db_dz_w1_zavg = mean(db_dz_w1(:,zzidx),2);
% % % % db_dz_w2_zavg = mean(db_dz_w2(:,zzidx),2);
% % % % db_dz_u_zavg = mean(db_dz_u(:,zzidx),2);
% % % % db_dz_3_zavg = mean(db_dz_3(:,zzidx),2);
% % % % db_dz_5_zavg = mean(db_dz_5(:,zzidx),2);
% % % % db_dz_reconst_zavg = mean(db_dz_reconst(:,zzidx),2);
% % % % dbdz_zavg = mean(dbdz(:,zzidx),2);
% % % % 
% % % % figure();clf;set(gcf,'color','w');
% % % % l6=plot(plot_time,db_dz_reconst_zavg(pidx),'LineWidth',1.5,'Color','k');
% % % % hold on;grid on;grid minor;box on;
% % % % l1=plot(plot_time,db_dz_w1_zavg(pidx),'LineWidth',1.5);
% % % % % l2=plot(plot_time,db_dz_w2_zavg(pidx),'LineWidth',1.5);
% % % % l3=plot(plot_time,db_dz_u_zavg(pidx),'LineWidth',1.5);
% % % % % l4=plot(plot_time,db_dz_3_zavg(pidx),'LineWidth',1.5);
% % % % % l5=plot(plot_time,db_dz_5_zavg(pidx),'LineWidth',1.5);
% % % % % plot(plot_time,dbdz_zavg(pidx),'LineWidth',1.5);
% % % % hold off;
% % % % 
% % % % % % % db_dz_leveln = dbdz(:,30)';
% % % % % % db_dz_leveln = db_dz_reconst(:,30)';
% % % % % % 
% % % % % % 
% % % % % % figure(4)
% % % % % % clf;set(gcf,'color','w');
% % % % % % hold on;grid on;grid minor;box on;
% % % % % % plot(plot_time,dB0_dz(pidx),'--','LineWidth',1.5);
% % % % % % plot(plot_time,dB_dz(pidx),'--','LineWidth',1.5);
% % % % % % plot(plot_time,db_dz_leveln(pidx),'--','LineWidth',1.5);
% % % % % % plot(plot_time,dB0_dz(pidx)+dB_dz(pidx)+db_dz_leveln(pidx),'LineWidth',2,'Color','k');
% % % % % % hold off;
% % % % % % set(gca,'FontSize',fontsize)
% % % % % % xlabel('Time (tidal cycles)');ylabel('N^2_{total} (s^{-2})')
% % % % % % xlim('tight')