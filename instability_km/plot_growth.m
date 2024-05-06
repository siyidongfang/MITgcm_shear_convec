
clear;
close all;

expdir = 'exps_flat/';
shear_all = [0:0.1e-4:1.8e-3];
rw_all= 10.^([-2:0.1:-1.1 -1.05 -1:0.0025:-0.7 -0.6:0.1:1]);
growth = zeros(length(shear_all),length(rw_all));
for i=1:length(shear_all)
    shear = shear_all(i);
    load([expdir '/growth_shear' num2str(shear*1e3,3) '_m01.mat'],'grow','topo','omega')
    growth(i,:)=grow;
    [max_growth(i) rw_idx]=max(grow);
    rw_mg(i) = rw_all(rw_idx);
end

criticalShear = omega*cosd(topo)/sind(topo);

plot_rw = atand(1./rw_all);
plot_rwmg = atand(1./rw_mg);

figure(22);
clf;
set(gcf,'Color','w');
pcolor(shear_all,plot_rw,growth')
shading flat;colormap(WhiteBlueGreenYellowRed(0))
hold on;
scatter(shear_all,plot_rwmg,10,'filled','black')
% plot(criticalShear*ones(1,nr),log10(rw_all),'Color','k','LineWidth',2)
grid on;grid minor;
title('Growth rate (1/hour)')
xlabel('Shear (1/s)')
ylabel('arctan(m_0/k) (degree)')
% ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
set(gca,'fontsize',20)
colorbar;
clim([0 6]/20)

figure(23)
clf;set(gcf,'Color','w')
plot(shear_all,max_growth,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
xlabel('Shear (1/s)')
title('Maximum growth rate (1/hour)')
ylabel('(1/hour)')

shear_km = shear_all;
growth_km = max_growth;
save('GrowthRate_km_topo4.mat','shear_km','growth_km','rw_mg')


figure(2)
clf;set(gcf,'Color','w')
plot(shear_all,1./rw_mg,'LineWidth',2)
xlabel('Shear (1/s)')
ylabel('m/k at t=0')
title('Aspect ratio of the initial perturbation (m_0/k) that corresponds to the largest growth rate')
set(gca,'fontsize',20)
grid on;grid minor;




% %%% Find out the wavenumber ratio rw=kx/mz corresponding to the maximum
% %%% growth rate, for each shear value
% growth_round = growth;
% growth_round(growth>1) = round(growth(growth>1),2);

%%% Plot timeseries of dbdz, u, w, ke, etc., at wavenumber ratio rw=kx/mz 
%%% corresponding to the maximum growth rate, for each shear value
%%% Analyze the spectra of w, dbdz, ke, etc.
%%
% for i=1:ns
for i=11
    showfig_dbdz = true;
    load(['output/topo' num2str(topo) '/shear' num2str(shear_all(i)) '_rw' num2str((rw_idx(i))) '.mat'])
    plot_timeseires;
end

%%
N2 = N^2;
if(plotRi)
    dt_ri = dt/1000;
    tt_ri = dt_ri:dt_ri:Nt*dt;
    isConvec = zeros(1,ns);
    Ri_min = zeros(1,ns);
    for i=1:ns 
        shear = shear_all(i);
        Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N2*cosd(topo) - N2*sind(topo)/omega*shear*sin(omega*tt_ri));
        % figure(50)
        % plot(tt_ri,Ri_inverse)
        if(min(Ri_inverse)<0)
            isConvec(i) = 1;
        end
        Ri_min(i) = 1/max(Ri_inverse);  
    end
    Ri_min(isConvec==1)=NaN;
    
    figure(40)
    plot(shear_all,(1./Ri_min))
    Ri_min(Ri_min==Inf)=NaN;
    

    figure(20);
    clf;
    set(gcf,'Color','w');
    pcolor(1./Ri_min,atand(1./rw_all),growth')
    shading flat;colormap(WhiteBlueGreenYellowRed(0))
    hold on;clim([0 6])
    grid on;grid minor;
    title('Growth rate (1/hour)')
    xlabel('1/Ri_{min}')
    ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
    set(gca,'fontsize',20)
    colorbar;
    
    
    figure(21)
    clf;set(gcf,'Color','w')
    plot(1./Ri_min,max_growth,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
    xlabel('1/Ri_{min}')
    title('Maximum growth rate (1/hour)')
    ylabel('(1/hour)')
end



