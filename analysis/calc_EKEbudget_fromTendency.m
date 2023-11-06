%%%
%%% calc_EKEbudget_fromTendency.m
%%%
%%% Calculate the Eddy Kinetic Energy budget from momentum tendency
%%% Since you're not using periodic buondary condition in the y-direction,
%%% remember to avoid the boundary point when calculating domain average or
%%% making plots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load the time-averaged data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prodir expname '_tavg_5days.mat'])
o1=1;
o2=24*5;

%%% Grid spacing matrices
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);

%%% load iteration numbers
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Time-mean momentum balance and velocities
uu = squeeze(UVEL);
vv = squeeze(VVEL);
ww = squeeze(WVEL);

Vm_Advec = squeeze(Vm_Advec);
Vm_dPhiY = squeeze(Vm_dPhiY);
Vm_Diss = squeeze(Vm_Diss);
Vm_Ext = squeeze(Vm_Ext);
AB_gV = squeeze(AB_gV);
Vm_ImplD = squeeze(Vm_ImplD);
tend_v = squeeze(TOTVTEND)/86400;

Um_Advec = squeeze(Um_Advec);
Um_dPhiX = squeeze(Um_dPhiX);
Um_Diss = squeeze(Um_Diss);
Um_Ext = squeeze(Um_Ext);
AB_gU = squeeze(AB_gU);
Um_ImplD = squeeze(Um_ImplD);
tend_u = squeeze(TOTUTEND)/86400;

Wm_Advec = squeeze(Wm_Advec);
Wm_Diss = squeeze(Wm_Diss);
AB_gW = squeeze(AB_gW);
tend_w = Wm_Advec+Wm_Diss+AB_gW; %%% not 100% sure

Vm_total = Vm_Advec+Vm_dPhiY+Vm_Diss+Vm_Ext+AB_gV+Vm_ImplD-tend_v;
Um_total = Um_Advec+Um_dPhiX+Um_Diss+Um_Ext+AB_gU+Um_ImplD-tend_u;

vv_advec_mean = vv.*Vm_Advec;
vv_dphi_mean = vv.*Vm_dPhiY;
vv_diss_mean = vv.*Vm_Diss;
vv_ext_mean = vv.*Vm_Ext;
vv_ab_mean = vv.*AB_gV;
vv_impld_mean = vv.*Vm_ImplD;
vv_tend_mean = vv.*tend_v;

uu_advec_mean = uu.*Um_Advec;
uu_dphi_mean = uu.*Um_dPhiX;
uu_diss_mean = uu.*Um_Diss;
uu_ext_mean = uu.*Um_Ext;
uu_ab_mean = uu.*AB_gU;
uu_impld_mean = uu.*Um_ImplD;
uu_tend_mean = uu.*tend_u;

ww_advec_mean = ww.*Wm_Advec;
ww_diss_mean = ww.*Wm_Diss;
ww_ab_mean = ww.*AB_gW;
ww_tend_mean = ww.*tend_w;

%%% Define matrices for estimating the total KE budget terms from
%%% high-frequency output
vv_advec = zeros(Nx,Nr);
vv_dphi = zeros(Nx,Nr);
vv_diss = zeros(Nx,Nr);
vv_ext = zeros(Nx,Nr);
vv_ab = zeros(Nx,Nr);
vv_impld = zeros(Nx,Nr);
vv_tend = zeros(Nx,Nr);

uu_advec = zeros(Nx,Nr);
uu_dphi = zeros(Nx,Nr);
uu_diss = zeros(Nx,Nr);
uu_ext = zeros(Nx,Nr);
uu_ab = zeros(Nx,Nr);
uu_impld = zeros(Nx,Nr);
uu_tend = zeros(Nx,Nr);

ww_advec = zeros(Nx,Nr);
ww_diss = zeros(Nx,Nr);
ww_ab = zeros(Nx,Nr);
ww_tend = zeros(Nx,Nr);


for o=o1:o2
    o
    nIter = dumpIters(o);

    vv_hour = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    uu_hour = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    ww_hour = squeeze(rdmds([exppath,'/results/WVEL'],nIter));

    Vm_Advec_hour = squeeze(rdmds([exppath,'/results/Vm_Advec'],nIter));
    Vm_dPhiY_hour = squeeze(rdmds([exppath,'/results/Vm_dPhiY'],nIter));
    Vm_Diss_hour = squeeze(rdmds([exppath,'/results/Vm_Diss'],nIter));
    Vm_Ext_hour = squeeze(rdmds([exppath,'/results/Vm_Ext'],nIter));
    AB_gV_hour = squeeze(rdmds([exppath,'/results/AB_gV'],nIter));
    Vm_ImplD_hour = squeeze(rdmds([exppath,'/results/Vm_ImplD'],nIter));
    tend_v_hour = squeeze(rdmds([exppath,'/results/TOTVTEND'],nIter))/86400;
    
    Um_Advec_hour = squeeze(rdmds([exppath,'/results/Um_Advec'],nIter));
    Um_dPhiX_hour = squeeze(rdmds([exppath,'/results/Um_dPhiX'],nIter));
    Um_Diss_hour = squeeze(rdmds([exppath,'/results/Um_Diss'],nIter));
    Um_Ext_hour = squeeze(rdmds([exppath,'/results/Um_Ext'],nIter));
    AB_gU_hour = squeeze(rdmds([exppath,'/results/AB_gU'],nIter));
    Um_ImplD_hour = squeeze(rdmds([exppath,'/results/Um_ImplD'],nIter));
    tend_u_hour = squeeze(rdmds([exppath,'/results/TOTUTEND'],nIter))/86400;
    
    Wm_Advec_hour = squeeze(rdmds([exppath,'/results/Wm_Advec'],nIter));
    Wm_Diss_hour = squeeze(rdmds([exppath,'/results/Wm_Diss'],nIter));
    AB_gW_hour = squeeze(rdmds([exppath,'/results/AB_gW'],nIter));
    tend_w_hour = Wm_Advec_hour+Wm_Diss_hour+AB_gW_hour; %%% not 100% sure

    vv_advec = vv_advec + vv_hour.*Vm_Advec_hour;
    vv_dphi = vv_dphi + vv_hour.*Vm_dPhiY_hour;
    vv_diss = vv_diss + vv_hour.*Vm_Diss_hour;
    vv_ext = vv_ext + vv_hour.*Vm_Ext_hour;
    vv_ab = vv_ab + vv_hour.*AB_gV_hour;
    vv_impld = vv_impld + vv_hour.*Vm_ImplD_hour;
    vv_tend = vv_tend + vv_hour.*tend_v_hour; 

    uu_advec = uu_advec + uu_hour.*Um_Advec_hour;
    uu_dphi = uu_dphi + uu_hour.*Um_dPhiX_hour;
    uu_diss = uu_diss + uu_hour.*Um_Diss_hour;
    uu_ext = uu_ext + uu_hour.*Um_Ext_hour;
    uu_ab = uu_ab + uu_hour.*AB_gU_hour;
    uu_impld = uu_impld + uu_hour.*Um_ImplD_hour;
    uu_tend = uu_tend + uu_hour.*tend_u_hour; 

    ww_advec = ww_advec + ww_hour.*Wm_Advec_hour;
    ww_diss = ww_diss + ww_hour.*Wm_Diss_hour;
    ww_ab = ww_ab + ww_hour.*AB_gW_hour;
    ww_tend = ww_tend + ww_hour.*tend_w_hour;

end

Nhours = o2-o1+1;

vv_advec = vv_advec/Nhours;
vv_dphi = vv_dphi/Nhours;
vv_diss = vv_diss/Nhours;
vv_ext = vv_ext/Nhours;
vv_ab = vv_ab/Nhours;
vv_impld = vv_impld/Nhours;
vv_tend = vv_tend/Nhours;

uu_advec = uu_advec/Nhours;
uu_dphi = uu_dphi/Nhours;
uu_diss = uu_diss/Nhours;
uu_ext = uu_ext/Nhours;
uu_ab = uu_ab/Nhours;
uu_impld = uu_impld/Nhours;
uu_tend = uu_tend/Nhours;

ww_advec = ww_advec/Nhours;
ww_diss = ww_diss/Nhours;
ww_ab = ww_ab/Nhours;
ww_tend = ww_tend/Nhours;

%%% Calculate the EKE budget
vv_advec_eddy = vv_advec_mean - vv_advec;
vv_dphi_eddy = vv_dphi_mean - vv_dphi;
vv_diss_eddy = vv_diss_mean - vv_diss;
vv_ext_eddy = vv_ext_mean - vv_ext;
vv_ab_eddy = vv_ab_mean - vv_ab;
vv_impld_eddy = vv_impld_mean - vv_impld;
vv_tend_eddy = vv_tend_mean - vv_tend;

uu_advec_eddy = uu_advec_mean - uu_advec;
uu_dphi_eddy = uu_dphi_mean - uu_dphi;
uu_diss_eddy = uu_diss_mean - uu_diss;
uu_ext_eddy = uu_ext_mean - uu_ext;
uu_ab_eddy = uu_ab_mean - uu_ab;
uu_impld_eddy = uu_impld_mean - uu_impld;
uu_tend_eddy = uu_tend_mean - uu_tend;

ww_advec_eddy = ww_advec_mean - ww_advec;
ww_diss_eddy = ww_diss_mean - ww_diss;
ww_ab_eddy = ww_ab_mean - ww_ab;
ww_tend_eddy = ww_tend_mean - ww_tend;


%%% from v grid to t grid
vv_advec_eddy_ugrid = (vv_advec_eddy+vv_advec_eddy([Ny 1:Ny-1],:))/2;
vv_dphi_eddy_ugrid = (vv_dphi_eddy+vv_dphi_eddy([Ny 1:Ny-1],:))/2;
vv_diss_eddy_ugrid = (vv_diss_eddy+vv_diss_eddy([Ny 1:Ny-1],:))/2;
vv_ext_eddy_ugrid = (vv_ext_eddy+vv_ext_eddy([Ny 1:Ny-1],:))/2;
vv_ab_eddy_ugrid = (vv_ab_eddy+vv_ab_eddy([Ny 1:Ny-1],:))/2;
vv_impld_eddy_ugrid = (vv_impld_eddy+vv_impld_eddy([Ny 1:Ny-1],:))/2;
vv_tend_eddy_ugrid = (vv_tend_eddy+vv_tend_eddy([Ny 1:Ny-1],:))/2;

eke_adv = vv_advec_eddy_ugrid + uu_advec_eddy;
eke_dphi = vv_dphi_eddy_ugrid + uu_dphi_eddy;
eke_diss = vv_diss_eddy_ugrid + uu_diss_eddy;
eke_ext = vv_ext_eddy_ugrid + uu_ext_eddy;
eke_ab = vv_ab_eddy_ugrid + uu_ab_eddy;
eke_impld = vv_impld_eddy_ugrid + uu_impld_eddy;
eke_tend = vv_tend_eddy_ugrid + uu_tend_eddy;

eke_adv = vv_advec_eddy_ugrid + uu_advec_eddy;
eke_dphi = vv_dphi_eddy_ugrid + uu_dphi_eddy;
eke_diss = vv_diss_eddy_ugrid + uu_diss_eddy;
eke_ext = vv_ext_eddy_ugrid + uu_ext_eddy;
eke_ab = vv_ab_eddy_ugrid + uu_ab_eddy;
eke_impld = vv_impld_eddy_ugrid + uu_impld_eddy;
eke_tend = vv_tend_eddy_ugrid + uu_tend_eddy;

eke_residual = eke_adv+eke_dphi+eke_diss+eke_ext+eke_ab+eke_impld-eke_tend;


eke_adv(eke_adv==0) = NaN;
eke_dphi(eke_dphi==0) = NaN;
eke_diss(eke_diss==0) = NaN;
eke_ext(eke_ext==0) = NaN;
eke_ab(eke_ab==0) = NaN;
eke_impld(eke_impld==0) = NaN;
eke_tend(eke_tend==0) = NaN;
eke_residual(eke_residual==0) = NaN;

save([prodir expname '_EKEbudget_fromTendency.mat'], ...
    'eke_adv','eke_dphi','eke_diss','eke_ext','eke_ab', ...
    'eke_impld','eke_tend','eke_residual')
%%
%%% Plot EKE budget

yy = xx;
[ZZ,YY] = meshgrid(zz,yy);
fontsize = 15;
CLIM = [-1 1]/1e7*3;
YLIM = [1000 1500];
figure(1)
clf;
set(gcf,'color','w','Position',[91 87 1038 993]);  
subplot(4,2,1)
pcolor(YY/1000,-ZZ,eke_adv); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim(YLIM);
title('EKE: Advection','interpreter','latex','FontSize',fontsize+4);

subplot(4,2,2)
pcolor(YY/1000,-ZZ,eke_dphi); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim(YLIM);
title('EKE: Pressure gradient force','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,3)
pcolor(YY/1000,-ZZ,eke_adv+eke_dphi); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim(YLIM);
title('EKE: Advec + Pres grad','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,4)
pcolor(YY/1000,-ZZ,eke_ext); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim(YLIM);
title('EKE: External forcing','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,5)
pcolor(YY/1000,-ZZ,eke_tend); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);
ylim(YLIM);
title('EKE: Tendency','interpreter','latex','FontSize',fontsize+4);



subplot(4,2,6)
pcolor(YY/1000,-ZZ,eke_diss+eke_impld); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/20);
ylim(YLIM);
% title('EKE: Dissipation (Explicit+Implicit)','interpreter','latex','FontSize',fontsize+4);
title('EKE: Dissipation','interpreter','latex','FontSize',fontsize+4);


% subplot(4,2,6)
% pcolor(YY/1000,-ZZ,eke_impld); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% % xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% handle=colorbar;set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);colormap redblue;
% clim(CLIM/40);
% ylim(YLIM);
% title('EKE: Dissipation (Implicit)','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,7)
pcolor(YY/1000,-ZZ,eke_ab); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/20);
ylim(YLIM);
title('EKE: Adams-Bashforth','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,8)
pcolor(YY/1000,-ZZ,eke_residual); shading interp;axis ij;set(gca,'color',gray);
hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/1e7);
ylim(YLIM);
title('EKE: Residual','interpreter','latex','FontSize',fontsize+4);

set(gcf, 'InvertHardcopy', 'off')
% print('-djpeg','-r150',[figdir 'EKE_fromTendency.jpeg']);




    
