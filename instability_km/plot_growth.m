
clear;
% close all;

% load(['output/growth_topo0.mat'])
load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1_Nsq1e-6_topo4/growth_topo4_NOdiff.mat')

%%%
%%% Calculate the Richardson Number
%%%
plotRi = true;
NaN_numbers = sum(isnan(growth),'all')
growth_largerw = growth(:,end);

%%% Find out the wavenumber ratio rw=kx/mz corresponding to the maximum
%%% growth rate, for each shear value
growth_round = growth;
growth_round(growth>1) = round(growth(growth>1),2);
[max_growth rw_idx] = max(growth_round,[],2);
% [max_growth rw_idx] = max(growth,[],2);
rw_mg = rw_all(rw_idx);

criticalShear = omega*cosd(topo)/sind(topo);


figure(22);
clf;
set(gcf,'Color','w');
pcolor(shear_all,atand(1./rw_all),growth')
shading flat;colormap(WhiteBlueGreenYellowRed(0))
hold on;
scatter(shear_all,atand(1./rw_mg),50,'filled','black')
plot(criticalShear*ones(1,nr),log10(rw_all),'Color','k','LineWidth',2)
grid on;grid minor;
title('Growth rate (1/hour)')
xlabel('Shear (1/s)')
ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
set(gca,'fontsize',20)
colorbar;
clim([0 6])

figure(23)
clf;set(gcf,'Color','w')
plot(shear_all,max_growth,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
xlabel('Shear (1/s)')
title('Maximum growth rate (1/hour)')
ylabel('(1/hour)')

figure(24)
clf;set(gcf,'Color','w')
plot(shear_all,growth_largerw,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
xlabel('Shear (1/s)')
title('Growth rate (1/hour), large k/m')
ylabel('(1/hour)')


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
    
    figure(25)
    clf;set(gcf,'Color','w')
    plot(1./Ri_min,growth_largerw,'LineWidth',2);grid on;grid minor;set(gca,'fontsize',20)
    xlabel('1/Ri_{min}')
    title('Growth rate (1/hour), large k/m')
    ylabel('(1/hour)')
    
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

