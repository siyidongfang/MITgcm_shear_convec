%%% 
%%% calc_GrowthRate.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations


clear;close all;
ne=1;
load_all

No = nDumps;
tidx = 1:No;
Nt = length(tidx);
Hshear = 250;
dz = delR(end);
Nshear = round(Hshear/dz);
zidx = Nr-Nshear:Nr-1;
% zidx = 1:Nr;
% Nshear = length(zidx)
botZ =-1500;
hab_shear = zz(zidx)-botZ;


uu_zavg = zeros(12,Nr);
tt_zavg = zeros(12,Nr);
for o = 1:12
    nIter = dumpIters(o);
    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    uu_zavg(o,:) = mean(uu,1);
    tt_zavg(o,:) = mean(tt,1);
end

wsq_xzavg = zeros(1,Nt);
upsq_xzavg = zeros(1,Nt);


for o = tidx
    nIter = dumpIters(o);
    time_h(o) = nIter.*deltaT./3600;
    time_tidal(o) = time_h(o)/12;

    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));

    wsq = ww.^2;
    remainder = rem(o,12);
    if(remainder==0)
        remainder = 12;
    end
    uu_zavg0 = repmat(uu_zavg(remainder,:),[Nx 1]);
    upsq = (uu-uu_zavg0).^2;

    wsq_xzavg(o) = mean(wsq(:,zidx),'all');
    upsq_xzavg(o) = mean(upsq(:,zidx),'all');

end

figure(1)
plot(time_h/12,log(wsq_xzavg/2)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(upsq_xzavg/2)/2,'LineWidth',2);


%%% Calculate the growth rate
fit_span = 12*12+1:22*12;
xxplot = time_h;
yyplot = log(wsq_xzavg/2+upsq_xzavg/2)/2;
[pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
grow = pp(1)
[y_fit,delta_fit] = polyval(pp,xxplot,S);



