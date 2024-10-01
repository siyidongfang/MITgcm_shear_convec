
clear;close all;

% load_colors;
fontsize = 16;
addpath ../instability/products/
load('fig4/Ri_flat.mat')

% load('fig4/topo4_km_kv2e-4.mat','max_grow','shear_all')
load('../instability_km/exps_varyingN/N1e-3output.mat')
[max_grow ~]=max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

load('GrowthRate_tanh_ZeroCenter_dz2.mat')
grow_Floquet_0Center = growth_Floquet;
for i=1:length(shear_Floquet)
    [a(i) b(i)] = min(abs(shear_Floquet(i)-shear_calc_Ri));
    Ri_Flo_0Center(i) = Ri_min(b(i));
end
clear shear_Floquet growth_Floquet

figure(10)
plot(1./Ri_Flo_0Center,grow_Floquet_0Center)
hold on;
plot(1./Ri_km,max_grow)
xlim([0 4])

%%
load('GrowthRate_tanh_BottomCenter_dz2.mat')
grow_Floquet_tanh = growth_Floquet;
for i=1:length(shear_Floquet)
    [a(i) b(i)] = min(abs(shear_Floquet(i)-shear_calc_Ri));
    Ri_Flo_tan(i) = Ri_min(b(i));
end
clear shear_Floquet growth_Floquet
plot(1./Ri_Flo_tan,grow_Floquet_tanh,'--')


load('GrowthRate_exps_linear_dz0.5.mat')
grow_Floquet_linear = growth_Floquet;
% load('GrowthRate_exps_linear.mat')
% grow_Floquet_linear = GrowthRate_KE;
for i=1:length(shear_Floquet)
    [a(i) b(i)] = min(abs(shear_Floquet(i)-shear_calc_Ri));
    Ri_Flo_linear(i) = Ri_min(b(i));
end
clear shear_Floquet growth_Floquet
plot(1./Ri_Flo_linear,grow_Floquet_linear,':')



grid on;
grid minor;
% figure(1)
% plot(shear_Floquet,growth_Floquet)


% load('GrowthRate_exps_tanh_ZeroCenter_dz2.mat')
% figure(1)
% plot(shear_Floquet,growth_Floquet)


% plot
%%% 8 panels 

%--- Flat bottom: growth rate as a function of inverse Ri and log10(lambda_x)

%--- Sloping bottom: growth rate as a function of inverse Ri and log10(lambda_x)

%--- Flat bottom: maximum growth rate as a function of inverse Ri (MITgcm vs theory)

%--- Sloping bottom: maximum growth rate as a function of inverse Ri (MITgcm vs theory)

%--- Flat bottom: buoyancy budget

%--- Sloping bottom: buoyancy budget

%--- Flat bottom: vorticity budget

%--- Sloping bottom: vorticity budget


% print('-dpng','-r200',['fig_supp/FigS_disc_growth.png']);




