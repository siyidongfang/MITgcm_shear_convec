%%%
%%% calc_EKEbudget.m
%%%
%%% Calculate the Eddy Kinetic Energy budget by yourself
%%% Since you're not using periodic buondary condition in the y-direction,
%%% remember to avoid the boundary point when calculating domain average or
%%% making plots.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load the time-averaged data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([prodir expname '_tavg_5days.mat'])
o1=1;
% o2=2*24*5-24;
o2 = 120;

vv = VVEL;
uu = UVEL;
ww = WVEL;
usq = UVELSQ;
vsq = VVELSQ;
wsq = WVELSQ;
tt = THETA;
wt = WVELTH;
uv = UV_VEL_Z; % Meridional Transport of Zonal Momentum (m^2/s^2), vorticity grid, middle level
uw = WU_VEL;   % Vertical Transport of Zonal Momentum, ugrid, lower level
vw = WV_VEL;   % Vertical Transport of Meridional Momentum, vgrid, lower level

vp = VVELPHI; % Mass-Weight Transp of Pressure Pot.(p/rho) Anomaly (m^3/s^3)
pp_anom = PHIHYD+PHI_NH;

% du_dt = TOTUTEND;
% dv_dt = TOTVTEND;

%%% Grid spacing matrices
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
DZC = repmat(reshape(-diff(zz),[1 1 Nr-1]),[Nx Ny 1]);

%%%%%%%%%%%%%%%%%%%%%
%%% Calculate EKE %%%
%%%%%%%%%%%%%%%%%%%%%
usq_eddy = usq-uu.^2;
vsq_eddy = vsq-vv.^2;
wsq_eddy = wsq-ww.^2;
EKE = 0.5 * ( 0.5 * (usq_eddy(1:Nx,:,:) + usq_eddy([2:Nx 1],:,:)) + 0.5 * (vsq_eddy(:,1:Ny,:) + vsq_eddy(:,[2:Ny 1],:)) ); %%% mass point
MKE = 0.5 * ( 0.5 * (uu(1:Nx,:,:).^2 + uu([2:Nx 1],:,:).^2) + 0.5 * (vv(:,1:Ny,:).^2 + vv(:,[2:Ny 1],:).^2) );

%%%%%%%%%%%%%%%%%%%%
%%% EKE tendency %%%
%%%%%%%%%%%%%%%%%%%%
%%% Assume to be zero??

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Boundary transport %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
usq_eddy_vgrid = usq_eddy;
usq_eddy_vgrid(2:Nx,:,:) = (usq_eddy_vgrid(1:Nx-1,:,:)+usq_eddy_vgrid(2:Nx,:,:))/2; %%% mass-grid
usq_eddy_vgrid(:,2:Ny,:) = (usq_eddy_vgrid(:,1:Ny-1,:)+usq_eddy_vgrid(:,2:Ny,:))/2; %%% v-grid

eke_vgrid = 0.5*(usq_eddy_vgrid + vsq_eddy);
eke_wgrid = zeros(Nx,Ny,Nr);
eke_wgrid(:,:,2:Nr) = 0.5 * (EKE(:,:,1:Nr-1) + EKE(:,:,2:Nr));

vm_eke = vv.*eke_vgrid;
wm_eke = ww.*eke_wgrid;
%%% Note that you have ignored all the \partial_x terms
dvm_eke_dy = (vm_eke(:,[2:Ny 1],:) - vm_eke(:,1:Ny,:)) ./ DY;
dwm_eke_dz = zeros(Nx,Ny,Nr);
dwm_eke_dz(:,:,2:Nr) = - diff(wm_eke,1,3) ./ DZC;

pp_anom_vgrid = pp_anom;
pp_anom_vgrid(:,2:Ny,:) = 0.5*(pp_anom(:,1:Ny-1,:)+pp_anom(:,2:Ny,:));
vp_eddy = vp-vv.*pp_anom_vgrid; % vgrid
dvp_eddy_dy = (vp_eddy(:,[2:Ny 1],:) - vp_eddy(:,1:Ny,:)) ./ DY;


%%%%%%%%%%%%%%%%%%%
%%% MKE --> EKE %%%
%%%%%%%%%%%%%%%%%%%
%%% Mean momentum fluxes
uv_mean = (0.5 * (uu(:,1:Ny,:) + uu(:,[Ny 1:Ny-1],:))) .* (0.5 * (vv(1:Nx,:,:) + vv([Nx 1:Nx-1],:,:))); %%% Mean zonal/meridional momentum flux 
uw_mean = zeros(Nx,Ny,Nr);
uw_mean(:,:,2:Nr) = (0.5 * (uu(:,:,1:Nr-1) + uu(:,:,2:Nr))) .* (0.5 * (ww(1:Nx,:,2:Nr) + ww([Nx 1:Nx-1],:,2:Nr))); 
uw_mean(:,:,1) = uu(:,:,1) .* (0.5 * (ww(1:Nx,:,1) + ww([Nx 1:Nx-1],:,1))); %%% Mean vertical flux of zonal momentum
vw_mean = zeros(Nx,Ny,Nr);
vw_mean(:,:,2:Nr) = (0.5 * (vv(:,:,1:Nr-1) + vv(:,:,2:Nr))) .* (0.5 * (ww(:,1:Ny,2:Nr) + ww(:,[Ny 1:Ny-1],2:Nr))); 
vw_mean(:,:,1) = vv(:,:,1) .* (0.5 * (ww(:,1:Ny,1) + ww(:,[Ny 1:Ny-1],1))); %%% Mean vertical flux of meridional momentum

%%% Eddy momentum fluxes
uv_eddy = uv - uv_mean; %%% Eddy zonal/meridional momentum flux, vorticity grid, cell center
uw_eddy = uw - uw_mean; %%% Eddy vertical flux of zonal momentum, ugrid, lower level
vw_eddy = vw - vw_mean; %%% Eddy vertical flux of meridional momentum, vgrid, lower level

%%% Mean derivatives
du_dx = (uu([2:Nx 1],:,:) - uu(1:Nx,:,:)) ./ DX;
du_dy = (uu(:,1:Ny,:) - uu(:,[Ny 1:Ny-1],:)) ./ DY;
du_dz = zeros(Nx,Ny,Nr);
du_dz(:,:,2:Nr) = - diff(uu,1,3) ./ DZC; %%% N.B. Does not account for partial cells
dv_dx = (vv(1:Nx,:,:) - vv([Nx 1:Nx-1],:,:)) ./ DX;
dv_dy = (vv(:,[2:Ny 1],:) - vv(:,1:Ny,:)) ./ DY;
dv_dz = zeros(Nx,Ny,Nr);
dv_dz(:,:,2:Nr) = - diff(vv,1,3) ./ DZC; %%% N.B. Does not account for partial cells

%%% Conversion terms
usq_du_dx = (0.5 * (usq_eddy(1:Nx,:,:) + usq_eddy([2:Nx 1],:,:))) .* du_dx; %%% mass point, cell center

uv_du_dy = uv_eddy .* du_dy;
uv_du_dy = 0.25 * (uv_du_dy(1:Nx,1:Ny,:) + uv_du_dy([2:Nx 1],1:Ny,:) + uv_du_dy(1:Nx,[2:Ny 1],:) + uv_du_dy([2:Nx 1],[2:Ny 1],:)); %%% from vorticity point to mass point

uw_du_dz = uw_eddy .* du_dz;
uw_du_dz(:,:,1:Nr-1) = 0.25 * (uw_du_dz(1:Nx,:,1:Nr-1) + uw_du_dz([2:Nx 1],:,1:Nr-1) + uw_du_dz(1:Nx,:,2:Nr) + uw_du_dz([2:Nx 1],:,2:Nr)); 
uw_du_dz(:,:,Nr) = 0.25 * (uw_du_dz(1:Nx,:,Nr) + uw_du_dz([2:Nx 1],:,Nr)); 

uv_dv_dx = uv_eddy .* dv_dx;
uv_dv_dyx= 0.25 * (uv_dv_dx(1:Nx,1:Ny,:) + uv_dv_dx([2:Nx 1],1:Ny,:) + uv_dv_dx(1:Nx,[2:Ny 1],:) + uv_dv_dx([2:Nx 1],[2:Ny 1],:));

vsq_dv_dy = (0.5 * (vsq_eddy(:,1:Ny,:) + vsq_eddy(:,[2:Ny 1],:))) .* dv_dy;

vw_dv_dz = vw_eddy .* dv_dz;
vw_dv_dz(:,:,1:Nr-1) = 0.25 * (vw_dv_dz(:,1:Ny,1:Nr-1) + vw_dv_dz(:,[2:Ny 1],1:Nr-1) + vw_dv_dz(:,1:Ny,2:Nr) + vw_dv_dz(:,[2:Ny 1],2:Nr)); 
vw_dv_dz(:,:,Nr) = 0.25 * (vw_dv_dz(:,1:Ny,Nr) + vw_dv_dz(:,[2:Ny 1],Nr)); 

%%% Overall energy conversion
MKE_EKE = -  (usq_du_dx + uv_du_dy + uw_du_dz + uv_dv_dx + vsq_dv_dy + vw_dv_dz);


%%%%%%%%%%%%%%%%%%%
%%% EPE --> EKE %%%
%%%%%%%%%%%%%%%%%%%
wt_mean = zeros(Nx,Ny,Nr+1);
wt_eddy = zeros(Nx,Ny,Nr+1);
wt_mean(:,:,2:Nr) = ww(:,:,2:Nr) .* (0.5*(tt(:,:,1:Nr-1)+tt(:,:,2:Nr))); %%% Mean vertical heat flux on vertical faces
wt_eddy(:,:,1:Nr) = wt - wt_mean(:,:,1:Nr); %%% Eddy vertical heat flux on vertical faces
wt_mean = 0.5 * (wt_mean(:,:,1:Nr) + wt_mean(:,:,2:Nr+1)); %%% Mean vertical heat flux at cell centers
wt_eddy = 0.5 * (wt_eddy(:,:,1:Nr) + wt_eddy(:,:,2:Nr+1)); %%% Eddy vertical heat flux at cell centers

EPE_EKE = tAlpha * gravity * wt_eddy; %%% EPE->EKE at cell centers




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimate the following EKE budget terms from hourly output 
%%% dwp_eddy_dz, dvp_eddy_dy
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

vp = zeros(Nx,Ny,Nr);
wp = zeros(Nx,Ny,Nr);
v_eke = zeros(Nx,Ny,Nr);
w_eke = zeros(Nx,Ny,Nr);
u_Um_Diss = zeros(Nx,Ny,Nr);
v_Vm_Diss = zeros(Nx,Ny,Nr);
u_Um_ImplD = zeros(Nx,Ny,Nr);
v_Vm_ImplD = zeros(Nx,Ny,Nr);
u_AB_gU = zeros(Nx,Ny,Nr);
v_AB_gV = zeros(Nx,Ny,Nr);

vert_mix_x = zeros(Nx,Ny,Nr);
vert_mix_y = zeros(Nx,Ny,Nr);

hori_diff_x = zeros(Nx,Ny,Nr);
hori_diff_y = zeros(Nx,Ny,Nr);

diff_x = zeros(Nx,Ny,Nr);
diff_y = zeros(Nx,Ny,Nr);


for o=o1:o2
    o
    nIter = dumpIters(o);
    vv_hour = rdmds([exppath,'/results/VVEL'],nIter);
    uu_hour = rdmds([exppath,'/results/UVEL'],nIter);
    ww_hour = rdmds([exppath,'/results/WVEL'],nIter);
    pp_hyd_hour = rdmds([exppath,'/results/PHIHYD'],nIter);
    pp_nh_hour = rdmds([exppath,'/results/PHI_NH'],nIter);
    pp_anom_hour = pp_hyd_hour + pp_nh_hour;
    vsq_hour = rdmds([exppath,'/results/VVELSQ'],nIter);
    usq_hour = rdmds([exppath,'/results/UVELSQ'],nIter);
    wsq_hour = rdmds([exppath,'/results/WVELSQ'],nIter);
    Um_Diss_hour = rdmds([exppath,'/results/Um_Diss'],nIter);
    Um_ImplD_hour = rdmds([exppath,'/results/Um_ImplD'],nIter);
    AB_gU_hour = rdmds([exppath,'/results/AB_gU'],nIter);
    Vm_Diss_hour = rdmds([exppath,'/results/Vm_Diss'],nIter);
    Vm_ImplD_hour = rdmds([exppath,'/results/Vm_ImplD'],nIter);
    AB_gV_hour = rdmds([exppath,'/results/AB_gV'],nIter);

    v4d_hour = rdmds([exppath,'/results/VA4DSMAG'],nIter);


    pp_anom_hour_wgrid = zeros(Nx,Ny,Nr);
    pp_anom_hour_wgrid(:,:,2:Nr) = 0.5 * (pp_anom_hour(:,:,1:Nr-1) + pp_anom_hour(:,:,2:Nr));
    wp = wp + ww_hour.*pp_anom_hour_wgrid; 

    pp_anom_hour_vgrid = pp_anom_hour;
    pp_anom_hour_vgrid(:,2:Ny,:) = 0.5*(pp_anom_hour(:,1:Ny-1,:)+pp_anom_hour(:,2:Ny,:));
    vp = vp + vv_hour.*pp_anom_hour_vgrid; 

    usq_eddy_hour = usq_hour-uu_hour.^2;
    vsq_eddy_hour = vsq_hour-vv_hour.^2;
    eke_hour = 0.5 * ( 0.5 * (usq_eddy_hour(1:Nx,:,:) + usq_eddy_hour([2:Nx 1],:,:))...
        + 0.5 * (vsq_eddy_hour(:,1:Ny,:) + vsq_eddy_hour(:,[2:Ny 1],:)) ); %%% mass point

    usq_eddy_hour_vgrid = usq_eddy_hour;
    usq_eddy_hour_vgrid(2:Nx,:,:) = (usq_eddy_hour_vgrid(1:Nx-1,:,:)+usq_eddy_hour_vgrid(2:Nx,:,:))/2; %%% mass-grid
    usq_eddy_hour_vgrid(:,2:Ny,:) = (usq_eddy_hour_vgrid(:,1:Ny-1,:)+usq_eddy_hour_vgrid(:,2:Ny,:))/2; %%% v-grid

    eke_hour_vgrid = 0.5*(usq_eddy_hour_vgrid + vsq_eddy_hour);
    eke_hour_wgrid = zeros(Nx,Ny,Nr);
    eke_hour_wgrid(:,:,2:Nr) = 0.5 * (eke_hour(:,:,1:Nr-1) + eke_hour(:,:,2:Nr));

    v_eke = v_eke + vv_hour .* eke_hour_vgrid;
    w_eke = w_eke + ww_hour .* eke_hour_wgrid;

    u_Um_Diss  = u_Um_Diss + uu_hour .* Um_Diss_hour;
    v_Vm_Diss  = v_Vm_Diss + vv_hour .* Vm_Diss_hour;
    u_Um_ImplD = u_Um_ImplD + uu_hour .* Um_ImplD_hour;
    v_Vm_ImplD = v_Vm_ImplD + vv_hour .* Vm_ImplD_hour;
    u_AB_gU    = u_AB_gU + uu_hour .* AB_gU_hour;
    v_AB_gV    = v_AB_gV + vv_hour .* AB_gV_hour;

    %%%% Vertical mixing
    dz2du_hour = zeros(Nx,Ny,Nr);
    dz2dv_hour = zeros(Nx,Ny,Nr);
    dz2du_hour(:,:,2:Nr) = - diff(uu_hour,1,3) ./ DZC; %%% N.B. Does not account for partial cells
    dz2du_hour(:,:,2:Nr-1) = - diff(dz2du_hour(:,:,2:Nr),1,3) ./ DZ(:,:,2:Nr-1); 
    dz2dv_hour(:,:,2:Nr) = - diff(vv_hour,1,3) ./ DZC; %%% N.B. Does not account for partial cells
    dz2dv_hour(:,:,2:Nr-1) = - diff(dz2dv_hour(:,:,2:Nr),1,3) ./ DZ(:,:,2:Nr-1); 
    vert_mix_x = vert_mix_x + uu_hour.*viscAr.*dz2du_hour; %%% u-grid, cell center
    vert_mix_y = vert_mix_y + vv_hour.*viscAr.*dz2dv_hour; %%% v-grid, cell center

    %%%% Horizontal diffusion
    d4y_uu_hour = zeros(Nx,Ny,Nr);
    d4y_vv_hour = zeros(Nx,Ny,Nr);

    d4y_uu_hour(:,3:Ny-2,:) = diff(uu_hour,4,2)./(delY(1).^4);
    d4y_vv_hour(:,3:Ny-2,:) = diff(vv_hour,4,2)./(delY(1).^4);

    d4y_uu_hour_tgrid = zeros(Nx,Ny,Nr);
    d4y_vv_hour_tgrid = zeros(Nx,Ny,Nr);
    d4y_uu_hour_tgrid(2:Nx,:,:) = (d4y_uu_hour(1:Nx-1,:,:)+ d4y_uu_hour(2:Nx,:,:))/2; %%% mass-grid
    d4y_vv_hour_tgrid(:,1:Ny-1,:) = (d4y_vv_hour(:,1:Ny-1,:)+d4y_vv_hour(:,2:Ny,:))/2;

    diff_x = diff_x + v4d_hour.*d4y_uu_hour_tgrid;
    diff_y = diff_y + v4d_hour.*d4y_vv_hour_tgrid;

    ud4yu_hour = uu_hour.*d4y_uu_hour;
    vd4yv_hour = vv_hour.*d4y_vv_hour;

    ud4yu_hour_tgrid = zeros(Nx,Ny,Nr);
    vd4yv_hour_tgrid = zeros(Nx,Ny,Nr);
    ud4yu_hour_tgrid(2:Nx,:,:) = (ud4yu_hour(1:Nx-1,:,:)+ ud4yu_hour(2:Nx,:,:))/2; %%% mass-grid
    vd4yv_hour_tgrid(:,1:Ny-1,:) = (vd4yv_hour(:,1:Ny-1,:)+vd4yv_hour(:,2:Ny,:))/2;

    hori_diff_x = hori_diff_x + v4d_hour.*ud4yu_hour_tgrid;
    hori_diff_y = hori_diff_y + v4d_hour.*vd4yv_hour_tgrid;


end

Nhours = o2-o1+1;
vp = vp/Nhours;
wp = wp/Nhours;
vert_mix_x = vert_mix_x/Nhours;
vert_mix_y = vert_mix_y/Nhours;
diff_x = diff_x/Nhours;
diff_y = diff_y/Nhours;
hori_diff_x = hori_diff_x/Nhours;
hori_diff_y = hori_diff_y/Nhours;
v_eke = v_eke/Nhours;
w_eke = w_eke/Nhours;
u_Um_Diss  = u_Um_Diss/Nhours;
v_Vm_Diss  = v_Vm_Diss/Nhours;
u_Um_ImplD = u_Um_ImplD/Nhours;
v_Vm_ImplD = v_Vm_ImplD/Nhours;
u_AB_gU    = u_AB_gU/Nhours;
v_AB_gV    = v_AB_gV/Nhours;


pp_anom_vgrid = pp_anom;
pp_anom_vgrid(:,2:Ny,:) = 0.5*(pp_anom(:,1:Ny-1,:)+pp_anom(:,2:Ny,:));
vp_eddy = vp-vv.*pp_anom_vgrid; % vgrid
dvp_eddy_dy = (vp_eddy(:,[2:Ny 1],:) - vp_eddy(:,1:Ny,:)) ./ DY;

pp_anom_wgrid = zeros(Nx,Ny,Nr);
pp_anom_wgrid(:,:,2:Nr) = 0.5 * (pp_anom(:,:,1:Nr-1) +  pp_anom(:,:,2:Nr));
wp_eddy = wp - ww.*pp_anom_wgrid;
dwp_eddy_dz = zeros(Nx,Ny,Nr);
dwp_eddy_dz(:,:,2:Nr) = - diff(wp_eddy,1,3) ./ DZC;

dv_eke_dy = (v_eke(:,[2:Ny 1],:) - v_eke(:,1:Ny,:)) ./ DY;
dw_eke_dz = zeros(Nx,Ny,Nr);
dw_eke_dz(:,:,2:Nr) = - diff(w_eke,1,3) ./ DZC;

bound_transp = - dvm_eke_dy - dwm_eke_dz   ...
               - dvp_eddy_dy - dwp_eddy_dz ...
               - dv_eke_dy - dw_eke_dz; % cell center, mass-grid


%%% dissipation terms???
u_Um_Diss_eddy = u_Um_Diss - uu.*Um_Diss;
v_Vm_Diss_eddy = v_Vm_Diss - vv.*Vm_Diss;
u_Um_ImplD_eddy = u_Um_ImplD - uu.*Um_ImplD;
v_Vm_ImplD_eddy = v_Vm_ImplD - vv.*Vm_ImplD;
u_AB_gU_eddy = u_AB_gU - uu.*AB_gU;
v_AB_gV_eddy = v_AB_gV - vv.*AB_gV;

%%%%%%%%%%%%%%%%%%%%%%%
%%% Vertical mixing %%%
%%%%%%%%%%%%%%%%%%%%%%%
dz2du = zeros(Nx,Ny,Nr);
dz2dv = zeros(Nx,Ny,Nr);
dz2du(:,:,2:Nr) = - diff(uu,1,3) ./ DZC; %%% N.B. Does not account for partial cells
dz2du(:,:,2:Nr-1) = - diff(dz2du(:,:,2:Nr),1,3) ./ DZ(:,:,2:Nr-1); 
dz2dv(:,:,2:Nr) = - diff(vv,1,3) ./ DZC; %%% N.B. Does not account for partial cells
dz2dv(:,:,2:Nr-1) = - diff(dz2dv(:,:,2:Nr),1,3) ./ DZ(:,:,2:Nr-1); 
vert_mix_x_mean = uu.*viscAr.*dz2du;
vert_mix_y_mean = vv.*viscAr.*dz2dv;
vert_mix_x_eddy = vert_mix_x-vert_mix_x_mean;
vert_mix_y_eddy = vert_mix_y-vert_mix_y_mean;

vert_mix_x_eddy_tgrid = zeros(Nx,Ny,Nr);
vert_mix_y_eddy_tgrid = zeros(Nx,Ny,Nr);
vert_mix_x_eddy_tgrid(2:Nx,:,:) = (vert_mix_x_eddy(1:Nx-1,:,:)+ vert_mix_x_eddy(2:Nx,:,:))/2; 
vert_mix_y_eddy_tgrid(:,1:Ny-1,:) = (vert_mix_y_eddy(:,1:Ny-1,:)+vert_mix_y_eddy(:,2:Ny,:))/2;

vert_mix_eddy = vert_mix_x_eddy_tgrid + vert_mix_y_eddy_tgrid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Horizontal diffusion %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu_tgrid = zeros(Nx,Ny,Nr);
vv_tgrid = zeros(Nx,Ny,Nr);
uu_tgrid(2:Nx,:,:) = (uu(1:Nx-1,:,:)+ uu(2:Nx,:,:))/2; 
vv_tgrid(:,1:Ny-1,:) = (vv(:,1:Ny-1,:)+vv(:,2:Ny,:))/2;

hori_diff_x_eddy = hori_diff_x-uu_tgrid.*diff_x;
hori_diff_y_eddy = hori_diff_y-vv_tgrid.*diff_y;

hori_diff_eddy = hori_diff_x_eddy + hori_diff_y_eddy;


v4z=squeeze(VA4ZSMAG); %%% SZ      MR|m^4/s           |Smagorinsky Biharm Visc Coeff. (m4/s) (Zeta Pt)
v4d=squeeze(VA4DSMAG); %%% SM      MR|m^4/s           |Smagorinsky Biharm Visc Coeff. (m4/s) (Div Pt)

v4z(v4z==0)=NaN;
v4d(v4d==0)=NaN;

%%
% yy=xx;
[ZZ,YY] = meshgrid(zz,yy);
figure(10)
clf
subplot(3,1,1)
pcolor(YY/1000,-ZZ/1000,v4z); shading interp;axis ij;colorbar;
clim([0 200])
subplot(3,1,2)
pcolor(YY/1000,-ZZ/1000,v4d); shading interp;axis ij;colorbar;
clim([0 200])
subplot(3,1,3)
pcolor(YY/1000,-ZZ/1000,v4z-v4d); shading interp;axis ij;colorbar;
colormap(redblue)
clim([-50 50])


EKE(EKE==0)=NaN;
EPE_EKE(EPE_EKE==0)=NaN;
MKE_EKE(MKE_EKE==0)=NaN;

EKE = squeeze(EKE);
EPE_EKE = squeeze(EPE_EKE);
MKE_EKE = squeeze(MKE_EKE);

bound_transp(bound_transp==0)=NaN; 
bound_transp = squeeze(bound_transp);

dvm_eke_dy(dvm_eke_dy==0)=NaN; 
dvm_eke_dy = squeeze(dvm_eke_dy);

dwm_eke_dz(dwm_eke_dz==0)=NaN; 
dwm_eke_dz = squeeze(dwm_eke_dz);

dvp_eddy_dy(dvp_eddy_dy==0)=NaN; 
dvp_eddy_dy = squeeze(dvp_eddy_dy);

dwp_eddy_dz(dwp_eddy_dz==0)=NaN; 
dwp_eddy_dz = squeeze(dwp_eddy_dz);

vert_mix_x(vert_mix_x==0)=NaN; 
vert_mix_x = squeeze(vert_mix_x);

vert_mix_y(vert_mix_y==0)=NaN; 
vert_mix_y = squeeze(vert_mix_y);

vert_mix_eddy(vert_mix_eddy==0)=NaN; 
vert_mix_eddy = squeeze(vert_mix_eddy);

hori_diff_x_eddy(hori_diff_x_eddy==0)=NaN; 
hori_diff_x_eddy = squeeze(hori_diff_x_eddy);

hori_diff_y_eddy(hori_diff_y_eddy==0)=NaN; 
hori_diff_y_eddy = squeeze(hori_diff_y_eddy);

hori_diff_eddy(hori_diff_eddy==0)=NaN; 
hori_diff_eddy = squeeze(hori_diff_eddy);


dv_eke_dy(dv_eke_dy==0)=NaN; 
dv_eke_dy = squeeze(dv_eke_dy);

dw_eke_dz(dw_eke_dz==0)=NaN; 
dw_eke_dz = squeeze(dw_eke_dz);

u_Um_Diss_eddy(u_Um_Diss_eddy==0)=NaN; 
u_Um_Diss_eddy = squeeze(u_Um_Diss_eddy);

v_Vm_Diss_eddy(v_Vm_Diss_eddy==0)=NaN; 
v_Vm_Diss_eddy = squeeze(v_Vm_Diss_eddy);

u_Um_ImplD_eddy(u_Um_ImplD_eddy==0)=NaN; 
u_Um_ImplD_eddy = squeeze(u_Um_ImplD_eddy);

v_Vm_ImplD_eddy(v_Vm_ImplD_eddy==0)=NaN; 
v_Vm_ImplD_eddy = squeeze(v_Vm_ImplD_eddy);

u_AB_gU_eddy(u_AB_gU_eddy==0)=NaN; 
u_AB_gU_eddy = squeeze(u_AB_gU_eddy);

v_AB_gV_eddy(v_AB_gV_eddy==0)=NaN; 
v_AB_gV_eddy = squeeze(v_AB_gV_eddy);




eke_dissipation = v_Vm_Diss_eddy + u_Um_Diss_eddy + v_Vm_ImplD_eddy + u_Um_ImplD_eddy + v_AB_gV_eddy+u_AB_gU_eddy;
eke_total = bound_transp + MKE_EKE + EPE_EKE + hori_diff_eddy + vert_mix_eddy + eke_dissipation;



%%%%%%%%%%%%%%%%%%%%%
%%% Save the data %%%
%%%%%%%%%%%%%%%%%%%%%

save([prodir expname '_EKEbudget_inst.mat'], ...
    'xx','yy','zz','YY','ZZ','EKE','eke_total',...
    'EPE_EKE','MKE_EKE','bound_transp','vert_mix_eddy','hori_diff_eddy','eke_dissipation',...
    'dvm_eke_dy','dwm_eke_dz','dvp_eddy_dy','dwp_eddy_dz','dv_eke_dy','dw_eke_dz', ...
    'hori_diff_x_eddy','hori_diff_y_eddy','vert_mix_x_eddy_tgrid','vert_mix_y_eddy_tgrid',...
    'v_Vm_Diss_eddy','u_Um_Diss_eddy','v_Vm_ImplD_eddy','u_Um_ImplD_eddy','v_AB_gV_eddy','u_AB_gU_eddy'...
    )


% %%
% 
% CLIM = [-1 1]*1e-3;
% 
% figure(2)
% clf;
% set(gcf,'color','w','Position',[55 282 975 695]);  
% subplot(3,2,1)
% pcolor(YY/1000,-ZZ/1000,bound_transp); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('boundary transport ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(3,2,2)
% pcolor(YY/1000,-ZZ/1000,EPE_EKE); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('$\mathrm{PE}\to\mathrm{EKE}$ ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(3,2,3)
% pcolor(YY/1000,-ZZ/1000,MKE_EKE);   shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('$\mathrm{MKE}\to\mathrm{EKE}$ ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(3,2,4)
% pcolor(YY/1000,-ZZ/1000,eke_dissipation); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('Dissipation ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(3,2,5)
% pcolor(YY/1000,-ZZ/1000,eke_total); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('Total EKE budget ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% set(gcf, 'InvertHardcopy', 'off')
% print('-djpeg','-r150',[figdir 'eke_fig3.jpeg']);
% 
% 
% 
% figure(1)
% clf;
% set(gcf,'color','w','Position',[55 1 1866 976]);  
% subplot(4,4,1)
% pcolor(YY/1000,-ZZ/1000,log10(EKE)); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap jet;
% clim([-6 -2]);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% title('$\log_{10}(\mathrm{EKE})$ ($\mathrm{m}^2$/$\mathrm{s}^2$)','interpreter','latex','FontSize',fontsize+2);
% 
% 
% subplot(4,4,2)
% pcolor(YY/1000,-ZZ/1000,EPE_EKE); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('$\mathrm{PE}\to\mathrm{EKE}$ ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(4,4,3)
% pcolor(YY/1000,-ZZ/1000,MKE_EKE);   shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('$\mathrm{MKE}\to\mathrm{EKE}$ ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(4,4,4)
% pcolor(YY/1000,-ZZ/1000,- dvm_eke_dy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('- dvm eke dy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(4,4,5)
% pcolor(YY/1000,-ZZ/1000,- dwm_eke_dz); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('- dwm eke dz ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,6)
% pcolor(YY/1000,-ZZ/1000,- dvp_eddy_dy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('- dvp eddy dy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(4,4,7)
% pcolor(YY/1000,-ZZ/1000,- dwp_eddy_dz); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('- dwp eddy dz ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(4,4,8)
% pcolor(YY/1000,-ZZ/1000,- dv_eke_dy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('dv eke dy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% 
% subplot(4,4,9)
% pcolor(YY/1000,-ZZ/1000, - dw_eke_dz); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('dw eke dz ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% % clim(CLIM);
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,10)
% pcolor(YY/1000,-ZZ/1000,bound_transp); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('boundary transport ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% 
% subplot(4,4,11)
% pcolor(YY/1000,-ZZ/1000,u_Um_Diss_eddy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('u UmDiss eddy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,12)
% pcolor(YY/1000,-ZZ/1000,v_Vm_Diss_eddy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('v VmDiss eddy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,13)
% pcolor(YY/1000,-ZZ/1000,u_Um_ImplD_eddy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('u UmImplD eddy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,14)
% pcolor(YY/1000,-ZZ/1000,v_Vm_ImplD_eddy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('v VmImplD eddy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,15)
% pcolor(YY/1000,-ZZ/1000,u_AB_gU_eddy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('u ABgU eddy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% subplot(4,4,16)
% pcolor(YY/1000,-ZZ/1000,v_AB_gV_eddy); shading interp;axis ij;set(gca,'color',gray);
% hold on;plot(yy/1000,(-bathy-40)/1000,'k--','LineWidth',1.5);plot(yy/1000,-bathy/1000,'k','LineWidth',1.5);hold off;
% xlabel('y (km)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (km)','interpreter','latex','FontSize',fontsize);
% title ('v ABgV eddy ($\mathrm{m}^2$/$\mathrm{s}^3$)','interpreter','latex','FontSize',fontsize+2)
% handle=colorbar;
% set(handle,'FontSize',fontsize);
% set(gca,'FontSize',fontsize);
% colormap redblue;
% clim(CLIM/1e4);
% ylim([1.5 -min(bathy)/1000+0.05])
% xlim([0.75 14.25])
% 
% set(gcf, 'InvertHardcopy', 'off')
% print('-djpeg','-r150',[figdir 'eke_fig4.jpeg']);

