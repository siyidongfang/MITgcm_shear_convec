%%%
%%% Compute velocity shear, given a Richardson number
%%%

clear


dt = 100;
NTtide = 100;
% omega = 1e-4;
omega = 2*pi/43200
Nt = NTtide/omega/dt;

topo=4;
N2 = 1e-6;

shear_convec = cosd(topo)/sind(topo)*omega;
if(topo==0)
    shear_convec = 5e-3;
end
shear_calc_Ri = 0:0.05e-5:(shear_convec+0.05e-5);
% shear_calc_Ri = 10e-4
% shear_calc_Ri = 0:0.25e-5:shear_convec;
ns = length(shear_calc_Ri);

%%
dt_ri = dt/1000;
tt_ri = dt_ri:dt_ri:Nt*dt;
isConvec = zeros(1,ns);
Ri_min = zeros(1,ns);
Ri_mean = zeros(1,ns);
Ri_harmonicmean = zeros(1,ns);
parfor i=1:ns 
    i
    shear = shear_calc_Ri(i);
    Ri_inverse = (cosd(topo)*shear*cos(omega*tt_ri)).^2./(N2 - N2*sind(topo)*cosd(topo)/omega*shear*sin(omega*tt_ri));
    % figure(50);set(gcf,'Position',[309 443 1200 400])
    % plot(tt_ri/43200*2*pi,Ri_inverse,'LineWidth',2);title('1/R_i(t)')
    % grid on;grid minor; set(gca,'FontSize',20);xlabel('Time (tidal cycles)')
    if(min(Ri_inverse)<0)
        isConvec(i) = 1;
    end
    Ri_min(i) = 1/max(Ri_inverse); 
    Ri_mean(i) = mean(1./Ri_inverse);
    % Ri_harmonicmean(i)=1/(sum(Ri_inverse)/length(Ri_inverse));
    Ri_harmonicmean(i) = harmmean(1./Ri_inverse);
end
% Ri_min(isConvec==1)=NaN;

% plot(shear_all,1./Ri_min)


[a idx] = min(abs(Ri_min-5))
Ri5 = Ri_min(idx);
shear_Ri5 = shear_calc_Ri(idx)


[a idx] = min(abs(Ri_min-1))
Ri1 = Ri_min(idx);
shear_Ri1 = shear_calc_Ri(idx)


[a idx] = min(abs(Ri_min-0.25));
Ri0_25 = Ri_min(idx);
shear_Ri0_25 = shear_calc_Ri(idx)

clear Ri_inverse tt_ri
% save('../figures/fig4/Ri_topo4_new.mat')
% save('../figures/fig4/Ri_flat_new.mat')

