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


for o=25*12+1:25*12+No

    nIter = dumpIters(o);
    time_h = nIter.*deltaT./3600;
    time_tidal(o) = time_h/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    tt(tt==0)=NaN;

    LC = 1*m1km;
    [mC,nC] = min(abs(xx-LC));
    ttC = tt(nC,:)+tt_background(nC,:);
    tt_timeseries(o,:) = tt(nC,:);
    rhoC = rhoConst.*(1-(ttC-tRef)*tAlpha);
    N2_timeseries(o,2:Nr) = -gravity/rhoConst.*(rhoC(1:end-1)-rhoC(2:end))./(zz(1:end-1)-zz(2:end));

    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    uu(uu==0)=NaN;
    ww(ww==0)=NaN;
    vv(vv==0)=NaN;
    S_timeseries(o,2:Nr) = (uu(nC,1:end-1)-uu(nC,2:end))./delR(2:end);
    uu_timeseries(o,:) = uu(nC,:);
    vv_timeseries(o,:) = vv(nC,:);
    ww_timeseries(o,:) = ww(nC,:);



end


botN = Nr;
botZ =zz(end);

% save('fig_supp/FigS_gcm_timeseries.mat','time_tidal','zz','botZ','tt_timeseries','N2_timeseries','uu_timeseries')

