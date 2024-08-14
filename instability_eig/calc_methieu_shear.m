%%% Plot growth rate as a function of omega0^2 and epsilon -- as the
%%% classical Arnold instability tongues


clear all;close all;

expdir = 'shear_methieu_noDiff/'
Ptide = 43200;
topo = 0;
N = 1e-3;
Ns = 1001;
shear_all = [0:3*N/(Ns-1):3*N]; 

for ns=1:Ns
    shear = shear_all(ns);
    load([expdir 'ptide' num2str(Ptide) '_topo' num2str(topo) '_N' num2str(N*1e3) '_shear' num2str(shear*1e3) '.mat']);
    grow(grow<=0)=NaN;
    grow = log(grow)/43200*3600;

    omega0 = N*kx_all/omega;

    figure(1)
    plot(lam_x_all,grow);
    xlim([-8000 8000])
    % plot(omega0,grow);
    % set(gca,'XScale','log')
    % xlim([0 100])
    ylim([-0.1 0.4])

    grow_ns(ns) = max(grow,[],'all','omitnan');
    grow_ns_positiverw(ns) = max(grow(8001:end),[],'all','omitnan');
    grow_ns_negativerw(ns) = max(grow(1:8000),[],'all','omitnan');
end



figure(3)
clf
plot(shear_all,grow_ns)
hold on;
plot(shear_all,grow_ns_positiverw,'--')
plot(shear_all,grow_ns_negativerw,':')

clear grow
save([expdir 'topo0_grow_ns.mat'])


