%%%
%%% calc_convec_probability.m
%%%
%%% Generate a map of N^<0 probability


clear; close all;
ne=1;
load_all;

Hz = sum(delR);
N2const = (1e-3)^2;
tNorth = N2const *(zz+Hz) /9.81/2e-4;
tt_background = ones(Nx,Nr);
for k=1:Nr
    tt_background(:,k) = squeeze(tt_background(:,k))*tNorth(k);
end

No = nDumps-1;
N2 = zeros(No,Nx,Nr);

for o=1:No
    nIter = dumpIters(o);
    tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
    tt = tt + tt_background;
    rho = rhoConst.*(1-(tt-tRef)*tAlpha);
    N2(o,:,1:Nr-1) = -gravity/rhoConst.*(rho(:,1:end-1)-rho(:,2:end))./(zz(1:end-1)-zz(2:end));
end

%%% Calculate N2<0 probability

convec_prob =  NaN*zeros(Nx,Nr);
for i=1:Nx
    for k=1:Nr-1
        convec_prob(i,k) = sum(N2(:,i,k)<0)/No;
    end
end

figure(1)
clf;set(gcf,'color','w');
pcolor(xx/1000,-zz,convec_prob');shading flat;colorbar;
axis ij;set(gca,'Fontsize',fontsize);
xlabel('x (km)','interpreter','latex');ylabel('Depth (m)','interpreter','latex');
title('Probability of $N^2<0$','interpreter','latex')
clim([0 0.3])
