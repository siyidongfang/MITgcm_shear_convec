%%%
%%% calc_EnergyBudget.m
%%%
%%% Convenience function to calculate energy budget terms.
%%% Use linear EOS instead of theta to calculate the buoyancy term.
%%%

    
clear; close all;

ne=1;
load_all;

%%% Load experiment data
load([prodir '/' expname '_tavg_5days.mat'],'UVEL','VVEL','WVEL','UVELSQ','VVELSQ','WVELSQ', ...
    'THETA','WVELTH','UV_VEL_Z','WU_VEL','WV_VEL'); %'SALT','WVELSLT'

vv=VVEL;
uu=UVEL;
ww=WVEL;
usq=UVELSQ;
vsq=VVELSQ;
wsq = WVELSQ;
tt=THETA;
wt=WVELTH;
uv = UV_VEL_Z; % Meridional Transport of Zonal Momentum (m^2/s^2)
uw = WU_VEL;   % Vertical Transport of Zonal Momentum
vw = WV_VEL;   % Vertical Transport of Meridional Momentum
tAlpha = 2e-4; % linear EOS thermal expansion coefficient (1/degC)

ss = zeros(Nx,Ny,Nr);
ws = zeros(Nx,Ny,Nr);
sBeta = 0;
zonalWind = zeros(Nx,Ny);
meridionalWind = zeros(Nx,Ny);
bottomDragLinear = 0;
% ss=SALT;
% ws=WVELSLT;
% sBeta = 7.4e-4; % linear EOS haline contraction coefficient (1/psu)
% load ([exppath '/input/setParams'],'Ua','Va')
% uwind = [Ua:-Ua/(Ny-1):0].*ones(Nx,1); 
% vwind = [Va:-Va/(Ny-1):0].*ones(Nx,1); 
% rho_a = 1.3;               %%% Air density, kg/m^3
% zonalWind = rho_a.*SEAICE_drag.*sqrt(uwind.^2+vwind.^2).*uwind;
% meridionalWind = rho_a.*SEAICE_drag.*sqrt(uwind.^2+vwind.^2).*vwind;


%%% Need the Eulerian-mean MOC to perform the buoyancy flux calculation
calcEulerianMOC;

%%% Grid spacing matrices
DX = repmat(delX',[1 Ny Nr]);
DY = repmat(delY,[Nx 1 Nr]);
DZ = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
DZC = repmat(reshape(-diff(zz),[1 1 Nr-1]),[Nx Ny 1]);

%%% Mesh grids for plotting
[YY,XX] = meshgrid(yy/1000,xx/1000);


%%%%%%%%%%%%%%%%%%%%%
%%% Calculate EKE %%%
%%%%%%%%%%%%%%%%%%%%%

usq_eddy = usq-uu.^2;
vsq_eddy = vsq-vv.^2;
wsq_eddy = wsq-ww.^2;
EKE = 0.5 * ( 0.5 * (usq_eddy(1:Nx,:,:) + usq_eddy([2:Nx 1],:,:)) + 0.5 * (vsq_eddy(:,1:Ny,:) + vsq_eddy(:,[2:Ny 1],:)) );
MKE = 0.5 * ( 0.5 * (uu(1:Nx,:,:).^2 + uu([2:Nx 1],:,:).^2) + 0.5 * (vv(:,1:Ny,:).^2 + vv(:,[2:Ny 1],:).^2) );
EKE_zavg = sum(EKE.*DZ.*hFacC,3) ./ sum(DZ.*hFacC,3); %%% Depth-averaged EKE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate PE->EKE,MKE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wt_mean = zeros(Nx,Ny,Nr+1);
wt_eddy = zeros(Nx,Ny,Nr+1);
wt_mean(:,:,2:Nr) = ww(:,:,2:Nr) .* (0.5*(tt(:,:,1:Nr-1)+tt(:,:,2:Nr))); %%% Mean vertical heat flux on vertical faces
wt_eddy(:,:,1:Nr) = wt - wt_mean(:,:,1:Nr); %%% Eddy vertical heat flux on vertical faces
wt_mean = 0.5 * (wt_mean(:,:,1:Nr) + wt_mean(:,:,2:Nr+1)); %%% Mean vertical heat flux at cell centers
wt_eddy = 0.5 * (wt_eddy(:,:,1:Nr) + wt_eddy(:,:,2:Nr+1)); %%% Eddy vertical heat flux at cell centers

ws_mean = zeros(Nx,Ny,Nr+1);
ws_eddy = zeros(Nx,Ny,Nr+1);
ws_mean(:,:,2:Nr) = ww(:,:,2:Nr) .* (0.5*(ss(:,:,1:Nr-1)+ss(:,:,2:Nr))); %%% Mean vertical heat flux on vertical faces
ws_eddy(:,:,1:Nr) = ws - ws_mean(:,:,1:Nr); %%% Eddy vertical heat flux on vertical faces
ws_mean = 0.5 * (ws_mean(:,:,1:Nr) + ws_mean(:,:,2:Nr+1)); %%% Mean vertical heat flux at cell centers
ws_eddy = 0.5 * (ws_eddy(:,:,1:Nr) + ws_eddy(:,:,2:Nr+1)); %%% Eddy vertical heat flux at cell centers

PE_MKE = tAlpha * gravity * wt_mean - sBeta * gravity * ws_mean; %%% PE->MKE at cell centers
PE_EKE = tAlpha * gravity * wt_eddy - sBeta * gravity * wt_eddy; %%% PE->EKE at cell centers

PE_EKE_zint = sum(PE_EKE.*DZ.*hFacC,3); %%% Depth-integrated PE->EKE
PE_EKE_zavg = PE_EKE_zint ./ sum(DZ.*hFacC,3); %%% Depth-averaged PE->EKE


%%% PE->MKE associated with zonal-mean velocity, temperature and salinity
w_xavg = ww;
w_xavg(w_xavg==0) = NaN;
w_xavg = squeeze(nanmean(w_xavg,1));
w_xavg(isnan(w_xavg)) = 0;
t_xavg = tt;
t_xavg(t_xavg==0) = NaN;
t_xavg = squeeze(nanmean(t_xavg,1));
t_xavg(isnan(t_xavg)) = 0;
s_xavg = ss;
s_xavg(s_xavg==0) = NaN;
s_xavg = squeeze(nanmean(s_xavg,1));
s_xavg(isnan(s_xavg)) = 0;
PE_ZKE_xyint = calcPE_ZKE_TS(hFacC,DX,DY,w_xavg,t_xavg,s_xavg,psiE,yy,gravity,tAlpha,sBeta);

%%% Integrate horizontally, starting from the southern edge of the
%%% Eulerian-mean streamfunction
PE_MKE_xyint = zeros(Nr,1);
PE_EKE_xyint = zeros(Nr,1);
for k=1:Nr
  psiE_k = 0.25*(psiE(1:Ny,k)+psiE(2:Ny+1,k)+psiE(1:Ny,k+1)+psiE(2:Ny+1,k+1));
  jrange = (yy>5e5) | (yy>4e5 & psiE_k'>0);
  PE_MKE_xyint(k) = squeeze(sum(sum(PE_MKE(:,jrange,k).*DX(:,jrange,k).*DY(:,jrange,k),1),2));
  PE_EKE_xyint(k) = squeeze(sum(sum(PE_EKE(:,jrange,k).*DX(:,jrange,k).*DY(:,jrange,k),1),2));
end

PE_MKE_xyint = PE_MKE_xyint - PE_ZKE_xyint;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate surface/bottom energy fluxes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dMKE_dt_wind = sum(sum(zonalWind(:,jrange).*uu(:,jrange,1).*DX(:,jrange,1).*DY(:,jrange,1)/rho0));
dEKE_dt_wind = 0;
isbot_mskW = ones(Nx,Ny,Nr);
isbot_mskW(:,:,1:Nr-1) = ceil(hFacW(:,:,1:Nr-1)) - ceil(hFacW(:,:,2:Nr));
isbot_mskW(:,:,Nr) = ceil(hFacW(:,:,Nr));
isbot_mskS = ones(Nx,Ny,Nr);
isbot_mskS(:,:,1:Nr-1) = ceil(hFacS(:,:,1:Nr-1)) - ceil(hFacS(:,:,2:Nr));
isbot_mskS(:,:,Nr) = ceil(hFacS(:,:,Nr));
dMKE_dt_bot = - sum(sum(sum(bottomDragLinear*uu(:,jrange,:).^2.*isbot_mskW(:,jrange,:).*DX(:,jrange,:).*DY(:,jrange,:)))) ...
              - sum(sum(sum(bottomDragLinear*vv(:,jrange,:).^2.*isbot_mskS(:,jrange,:).*DX(:,jrange,:).*DY(:,jrange,:))));
dEKE_dt_bot = - sum(sum(sum(bottomDragLinear*usq_eddy(:,jrange,:).*isbot_mskW(:,jrange,:).*DX(:,jrange,:).*DY(:,jrange,:)))) ...
              - sum(sum(sum(bottomDragLinear*vsq_eddy(:,jrange,:).*isbot_mskS(:,jrange,:).*DX(:,jrange,:).*DY(:,jrange,:))));


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate MKE->EKE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Mean momentum fluxes
uv_mean = (0.5 * (uu(:,1:Ny,:) + uu(:,[Ny 1:Ny-1],:))) .* (0.5 * (vv(1:Nx,:,:) + vv([Nx 1:Nx-1],:,:))); %%% Mean zonal/meridional momentum flux 
uw_mean = zeros(Nx,Ny,Nr);
uw_mean(:,:,2:Nr) = (0.5 * (uu(:,:,1:Nr-1) + uu(:,:,2:Nr))) .* (0.5 * (ww(1:Nx,:,2:Nr) + ww([Nx 1:Nx-1],:,2:Nr))); 
uw_mean(:,:,1) = uu(:,:,1) .* (0.5 * (ww(1:Nx,:,1) + ww([Nx 1:Nx-1],:,1))); %%% Mean vertical flux of zonal momentum
vw_mean = zeros(Nx,Ny,Nr);
vw_mean(:,:,2:Nr) = (0.5 * (vv(:,:,1:Nr-1) + vv(:,:,2:Nr))) .* (0.5 * (ww(:,1:Ny,2:Nr) + ww(:,[Ny 1:Ny-1],2:Nr))); 
vw_mean(:,:,1) = vv(:,:,1) .* (0.5 * (ww(:,1:Ny,1) + ww(:,[Ny 1:Ny-1],1))); %%% Mean vertical flux of meridional momentum

%%% Eddy momentum fluxes
uv_eddy = uv - uv_mean; %%% Eddy zonal/meridional momentum flux 
uw_eddy = uw - uw_mean; %%% Eddy vertical flux of zonal momentum
vw_eddy = vw - vw_mean; %%% Eddy vertical flux of meridional momentum

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
usq_du_dx = (0.5 * (usq_eddy(1:Nx,:,:) + usq_eddy([2:Nx 1],:,:))) .* du_dx;

uv_du_dy = uv_eddy .* du_dy;
uv_du_dy = 0.25 * (uv_du_dy(1:Nx,1:Ny,:) + uv_du_dy([2:Nx 1],1:Ny,:) + uv_du_dy(1:Nx,[2:Ny 1],:) + uv_du_dy([2:Nx 1],[2:Ny 1],:));

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
MKE_EKE_hori = -  (usq_du_dx + uv_du_dy + uv_dv_dx + vsq_dv_dy); %%% Horizontal component
MKE_EKE_vert = -  (uw_du_dz + vw_dv_dz);


EKE(EKE==0)=NaN;
PE_EKE(PE_EKE==0)=NaN;
MKE_EKE(MKE_EKE==0)=NaN;
PE_MKE(PE_MKE==0)=NaN;

EKE = squeeze(EKE);
PE_EKE = squeeze(PE_EKE);
MKE_EKE = squeeze(MKE_EKE);
PE_MKE = squeeze(PE_MKE);

%%% Store computed data for later
save([prodir '/' expname,'_EnergyBudget_3days.mat'],'DX','DY','DZ','DZC','YY','XX',...
    'EKE','PE_EKE','MKE_EKE','PE_MKE',...
    'MKE_EKE_hori','MKE_EKE_vert'); 



