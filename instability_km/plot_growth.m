

% load("growth_topo0.mat")

%%%
%%% Calculate the Richardson Number
%%%
N2 = N^2;
for i=1:ns 
    shear = shear_all(i);
    Ri_inverse = (shear*cos(tt)).^2./(N2*cosd(topo) - N2*sind(topo)/omega*shear*sin(tt));
    if(min(Ri_inverse)<0)
        isConvec(i) = 1;
    end
    Ri_min(i) = 1/max(Ri_inverse);
end

Ri_min(Ri_min==Inf)=NaN;
max_growth = max(growth,[],2)

figure(20);
clf;
set(gcf,'Color','w');
pcolor(1./Ri_min,log10(rw_all),growth')
shading flat;colormap(WhiteBlueGreenYellowRed(0))
hold on;
grid on;grid minor;
title('Growth rate (1/hour)')
xlabel('1/Ri_{min}')
ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
set(gca,'fontsize',20)
colorbar;

figure(21)
plot(1./Ri_min,max_growth)
xlabel('1/Ri_{min}')
title('Maximum growth rate (1/hour)')
ylabel('(1/hour)')

figure(22);
clf;
set(gcf,'Color','w');
pcolor(shear_all,log10(rw_all),growth')
shading flat;colormap(WhiteBlueGreenYellowRed(0))
hold on;
grid on;grid minor;
title('Growth rate (1/hour)')
xlabel('Shear (1/s)')
ylabel('Wavenumber ratio log_{10}(k_x/m_z)')
set(gca,'fontsize',20)
colorbar;

figure(23)
plot(shear_all,max_growth)
xlabel('Shear (1/s)')
title('Maximum growth rate (1/hour)')
ylabel('(1/hour)')
