addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
expdir = '../exps_topo4/';
loadexp;

rhoConst = 999.8;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);


No = 4*12;
tidx = 25*12+1:25*12+No;

uu_timeseries= zeros(No,Nr);
S_timeseries = zeros(No,Nr);
vv_timeseries = zeros(No,Nr);
ww_timeseries = zeros(No,Nr);
tt_timeseries = zeros(No,Nr);
N2_timeseries = zeros(No,Nr);
epsilon_timeseries = zeros(No,Nr);
chi_timeseries = zeros(No,Nr);

time_tidal = zeros(1,No);
pp_mid = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined

tRef = 0;
m1km = 1000;

Hz = sum(delR);
N2const = 1e-6;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);

for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end



for i=1:length(tidx)

    o=tidx(i);

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    time_tidal(o) = time_h/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    tt(tt==0)=NaN;

    LC = 1*m1km;
    [mC,nC] = min(abs(xx-LC));
    ttC = tt(nC,:)+tt_background(nC,:);
    tt_timeseries(i,:) = ttC;
    rhoC = rhoConst.*(1-(ttC-tRef)*tAlpha);
    N2_timeseries(i,2:Nr) = -gravity/rhoConst.*(rhoC(1:end-1)-rhoC(2:end))./(zz(1:end-1)-zz(2:end));

    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    uu(uu==0)=NaN;
    ww(ww==0)=NaN;
    vv(vv==0)=NaN;
    S_timeseries(i,2:Nr) = (uu(nC,1:end-1)-uu(nC,2:end))./delR(2:end);
    uu_timeseries(i,:) = uu(nC,:);
    vv_timeseries(i,:) = vv(nC,:);
    ww_timeseries(i,:) = ww(nC,:);


    dz = delR(1);
    dx = delX(1);
    du_dz = zeros(Nx,Nr);
    dv_dz = zeros(Nx,Nr);
    dw_dz = zeros(Nx,Nr);
    du_dx = zeros(Nx,Nr);
    dv_dx = zeros(Nx,Nr);
    dw_dx = zeros(Nx,Nr);
    w_tlevel = zeros(Nx,Nr);
    du_dx_ugrid_wlev = zeros(Nx,Nr);
    dv_dx_wlev = zeros(Nx,Nr);
    
    du_dz(:,2:Nr) = - diff(uu,1,2) ./dz; 
    dv_dz(:,2:Nr) = - diff(vv,1,2) ./dz; % on v-grid
    dv_dz_ugrid = 0.5.*(dv_dz([Nx 1:Nx-1],:,:)+dv_dz(1:Nx,:,:));
    w_tlevel(:,1:Nr-1) = 0.5*(ww(:,1:Nr-1)+ww(:,2:Nr));
    dw_dz(:,2:Nr) = - diff(w_tlevel,1,2) ./dz; % on v-grid
    dw_dz_ugrid = 0.5.*(dw_dz([Nx 1:Nx-1],:,:)+dw_dz(1:Nx,:,:));
    
    du_dx= (uu([2:Nx 1],:)-uu) ./dx; % on v-grid
    du_dx_ugrid = 0.5.*(du_dx([Nx 1:Nx-1],:,:)+du_dx(1:Nx,:,:));
    du_dx_ugrid_wlev(:,2:Nr) = 0.5*(du_dx_ugrid(:,1:Nr-1)+du_dx_ugrid(:,2:Nr));
    
    dw_dx= (ww-ww([Nx 1:Nx-1],:)) ./dx; 
    
    dv_dx= (vv-vv([Nx 1:Nx-1],:)) ./dx; 
    dv_dx_wlev(:,2:Nr)= 0.5*(dv_dx(:,1:Nr-1)+dv_dx(:,2:Nr));
    
    
    epsilon = viscAr.*(du_dz.^2+dv_dz_ugrid.^2+dw_dz_ugrid.^2  ...
        +du_dx_ugrid_wlev.^2+dv_dx_wlev.^2+dw_dx.^2);  %%% Mass-point, lower level
    
    dT_dx= (tt-tt([Nx 1:Nx-1],:)) ./dx; % u-grid, t-lev
    
    tt_ugrid = 0.5*(tt([Nx 1:Nx-1],:)+tt);
    tt_ugrid_wlev = zeros(Nx,Nr);
    tt_ugrid_wlev(:,2:Nr) =  0.5*(tt_ugrid(:,1:Nr-1)+tt_ugrid(:,2:Nr));
    
    dT_dz = zeros(Nx,Nr);
    dT_dz(:,1:Nr-1) = - diff(tt_ugrid_wlev,1,2) ./ dz;
    
    chi_x = diffKhT.*(dT_dx.^2); %%% u-grid, t-lev
    chi_z = diffKrT.*dT_dz.^2;   %%% u-grid, t-lev
    chi_tt = chi_z + chi_x; %%% u-grid, t-lev


    epsilon_timeseries(i,:) = epsilon(nC,:);
    chi_timeseries(i,:) = chi_tt(nC,:);

end

botN = Nr;
botZ =zz(end);

save('fig_supp/FigS_gcm_timeseries.mat','tidx',...
    'time_tidal','zz','botZ','tt_timeseries','N2_timeseries','uu_timeseries',...
    'vv_timeseries','ww_timeseries','S_timeseries','epsilon_timeseries','chi_timeseries')

