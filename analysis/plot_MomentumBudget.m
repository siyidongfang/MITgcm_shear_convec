%%%
%%% plot_MomentumBudget.m
%%%
%%% Convenience script to plot the momentum budget from momentum tendency diagnostics.
%%%
%%% TOTUTEND/86400=-g*d(Eta)/dx+Um_Diss+Um_Advec+Um_dPHdx+Um_Ext+AB_gU+VISrI_Um/volume
%%% The new diagnostics Um_dPhiX and Vm_dPhiY have included the terms -g*d(Eta)/dx and -g*d(Eta)/dy?
%%% Be sure to get the tendency from the vertical viscous flux with the following equation for the kth level:
%%% ( VISrI_Um^k+1 - VISrI_Um^k ) / Cell_Vol
%%% For the surface layer, VISrI_Um^1 should be zero.

% % % With implicit viscosity, you need the extra terms that are derived from VISrI_Um and Vm (the vertical diffusive fluxes of momentum, that are computed with an implicit method).
% % % 2) Um_Diss is the tendency due to dissipation, but VISrE/I_Um is the vertical flux (du/dz*area). To get to a tendency, you need to compute (VISrE_Um[k]-VISrE_Um[k+1])/volumeOfGridCell. 
% % % 3) the same is true for 'DFrE_TH' and ‘DFrI_TH’. They are vertical fluxes that already include the area. The units are correct to, as far as I can see, because the diffusion term (d/dz kappa dtheta/dz) should have degC/s, so kappa dtheta/dz has the unis degC/s * m, and kappa dtheta/dz * area then degC/s m^3, you need to do the same computation as for VISrE/I_Um to get the tendencies.


clear;close all;
ne=1;
load_all;
rho0 = rhoConst;
% load([prodir expname '_tavg_5days.mat'])

o1 = 1;
o2 = 48;
figdir = [expdir expname '/'];

for o=o1:o2
nIter = dumpIters(o);
time_h = nIter.*deltaT./3600;

Vm_Advec = rdmds([exppath,'/results/Vm_Advec'],nIter);
Vm_dPhiY = rdmds([exppath,'/results/Vm_dPhiY'],nIter);
Vm_Diss = rdmds([exppath,'/results/Vm_Diss'],nIter);
Vm_Ext = rdmds([exppath,'/results/Vm_Ext'],nIter);
AB_gV = rdmds([exppath,'/results/AB_gV'],nIter);
Vm_ImplD = rdmds([exppath,'/results/Vm_ImplD'],nIter);
TOTVTEND = rdmds([exppath,'/results/TOTVTEND'],nIter);

Um_Advec = rdmds([exppath,'/results/Um_Advec'],nIter);
Um_dPhiX = rdmds([exppath,'/results/Um_dPhiX'],nIter);
Um_Diss = rdmds([exppath,'/results/Um_Diss'],nIter);
Um_Ext = rdmds([exppath,'/results/Um_Ext'],nIter);
AB_gU = rdmds([exppath,'/results/AB_gU'],nIter);
Um_ImplD = rdmds([exppath,'/results/Um_ImplD'],nIter);
TOTUTEND = rdmds([exppath,'/results/TOTUTEND'],nIter);

Wm_Advec = rdmds([exppath,'/results/Wm_Advec'],nIter);
Wm_Diss = rdmds([exppath,'/results/Wm_Diss'],nIter);
AB_gW = rdmds([exppath,'/results/AB_gW'],nIter);

%%% Grid spacing matrices
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
DZC = repmat(reshape(-diff(zz),[1 1 Nr-1]),[Nx Ny 1]);

tend_v = TOTVTEND/86400;
tend_u = TOTUTEND/86400;

%%% The results from the calculation below exactly equal to Vm_ImplD and Um_ImplD 
% Cell_Vol = delX(1).*delY(1).*DZC;
% Vm_Visc = zeros(Nx,Ny,Nr); %%% tendency from the vertical viscous flux 
% Vm_Visc(:,:,2:end) = diff(VISrI_Vm,1,3)./Cell_Vol;
% Um_Visc = zeros(Nx,Ny,Nr); %%% tendency from the vertical viscous flux 
% Um_Visc(:,:,2:end) = diff(VISrI_Um,1,3)./Cell_Vol;

Vm_total = Vm_Advec+Vm_dPhiY+Vm_Diss+Vm_Ext+AB_gV+Vm_ImplD-tend_v;
Um_total = Um_Advec+Um_dPhiX+Um_Diss+Um_Ext+AB_gU+Um_ImplD-tend_u;
Wm_total = Wm_Advec+Wm_Diss+AB_gW;

Vm_total = squeeze(Vm_total);
Vm_Advec = squeeze(Vm_Advec);
Vm_dPhiY = squeeze(Vm_dPhiY);
Vm_Diss = squeeze(Vm_Diss);
Vm_Ext = squeeze(Vm_Ext);
% Vm_Visc = squeeze(Vm_Visc);
AB_gV = squeeze(AB_gV);
Vm_ImplD = squeeze(Vm_ImplD);
tend_v = squeeze(tend_v);

Um_total = squeeze(Um_total);
Um_Advec = squeeze(Um_Advec);
Um_dPhiX = squeeze(Um_dPhiX);
Um_Diss = squeeze(Um_Diss);
% Um_Visc = squeeze(Um_Visc);
AB_gU = squeeze(AB_gU);
Um_ImplD = squeeze(Um_ImplD);
tend_u = squeeze(tend_u);

Wm_total = squeeze(Wm_total);
Wm_Advec = squeeze(Wm_Advec);
Wm_Diss = squeeze(Wm_Diss);
AB_gW = squeeze(AB_gW);


if(Ny>1)
    [ZZ,YY] = meshgrid(zz,yy);
end
if(Nx>1)
    [ZZ,YY] = meshgrid(zz,xx);
    yy = xx;
end
fontsize = 15;
CLIM = [-2 2]/1e5;



%%
figure(1)
clf;
set(gcf,'color','w','Position',[91 87 1038 993]);  
subplot(4,2,1)
pcolor(YY/1000,-ZZ/1000,Vm_Advec); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Advection','interpreter','latex','FontSize',fontsize+4);

subplot(4,2,2)
pcolor(YY/1000,-ZZ/1000,Vm_dPhiY); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/1e20/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Pressure gradient force','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,3)
pcolor(YY/1000,-ZZ/1000,Vm_Advec+Vm_dPhiY); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Advec + Pres grad','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,4)
pcolor(YY/1000,-ZZ/1000,Vm_Diss); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/100000/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Dissipation (Explicit)','interpreter','latex','FontSize',fontsize+4);

subplot(4,2,5)
pcolor(YY/1000,-ZZ/1000,tend_v); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Tendency','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,6)
pcolor(YY/1000,-ZZ/1000,Vm_ImplD); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/100000);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Dissipation (Implicit)','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,7)
pcolor(YY/1000,-ZZ/1000,AB_gV); shading interp;axis ij;set(gca,'color',gray);
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/5000);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Adams-Bashforth','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,8)
pcolor(YY/1000,-ZZ/1000,Vm_total); shading interp;axis ij;set(gca,'color',gray);
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/2/1e7/1e7);ylim([1.1 -min(bathy)/1000]);
title('Along-isobath V: Residual','interpreter','latex','FontSize',fontsize+4);

set(gcf, 'InvertHardcopy', 'off')
print('-djpeg','-r150',[figdir 'evo2_V_momentum' num2str(o) '.jpeg']);

%%

% CLIM = CLIM/10;

figure(2)
clf;
set(gcf,'color','w','Position',[91 87 1038 993]);  

subplot(4,2,1)
pcolor(YY/1000,-ZZ/1000,Um_Advec); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Advection','interpreter','latex','FontSize',fontsize+4);

subplot(4,2,2)
pcolor(YY/1000,-ZZ/1000,Um_dPhiX); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/5e3/10/300);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Pressure gradient force','interpreter','latex','FontSize',fontsize+4);



subplot(4,2,3)
pcolor(YY/1000,-ZZ/1000,Um_Advec); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Advection','interpreter','latex','FontSize',fontsize+4);

subplot(4,2,4)
pcolor(YY/1000,-ZZ/1000,Um_Diss); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/1e5/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Dissipation (Explicit)','interpreter','latex','FontSize',fontsize+4);

subplot(4,2,5)
pcolor(YY/1000,-ZZ/1000,tend_u); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Tendency','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,6)
pcolor(YY/1000,-ZZ/1000,Um_ImplD); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/1e4);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Dissipation (Implicit)','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,7)
pcolor(YY/1000,-ZZ/1000,AB_gU); shading interp;axis ij;set(gca,'color',gray);
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/1000);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Adams-Bashforth','interpreter','latex','FontSize',fontsize+4);


subplot(4,2,8)
pcolor(YY/1000,-ZZ/1000,Um_total); shading interp;axis ij;set(gca,'color',gray);
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/2/1e7/4e6);ylim([1.1 -min(bathy)/1000]);
title('Cross-isobath U: Residual','interpreter','latex','FontSize',fontsize+4);


set(gcf, 'InvertHardcopy', 'off')
print('-djpeg','-r150',[figdir 'evo2_U_momentum' num2str(o) '.jpeg']);


figure(3)
clf;
set(gcf,'color','w','Position', [70 411 1058 477]);  
subplot(2,2,1)
pcolor(YY/1000,-ZZ/1000,Wm_Advec); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/500/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Vertical: Advection','interpreter','latex','FontSize',fontsize+2);

subplot(2,2,2)
pcolor(YY/1000,-ZZ/1000,Wm_Diss); shading interp;axis ij;set(gca,'color',gray);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/100000/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Vertical: Dissipation','interpreter','latex','FontSize',fontsize+2);

subplot(2,2,3)
pcolor(YY/1000,-ZZ/1000,AB_gW); shading interp;axis ij;set(gca,'color',gray);
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/1e4/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Vertical: Adams-Bashforth','interpreter','latex','FontSize',fontsize+4);

subplot(2,2,4)
pcolor(YY/1000,-ZZ/1000,Wm_total); shading interp;axis ij;set(gca,'color',gray);
xlabel('x (km)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
handle=colorbar;set(handle,'FontSize',fontsize);
set(gca,'FontSize',fontsize);colormap redblue;
clim(CLIM/500/10/1e4);ylim([1.1 -min(bathy)/1000]);
title('Vertical: Residual (==tendency?)','interpreter','latex','FontSize',fontsize+2);

set(gcf, 'InvertHardcopy', 'off')
print('-djpeg','-r150',[figdir 'evo2_W_momentum' num2str(o) '.jpeg']);


end
