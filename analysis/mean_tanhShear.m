
clear all;close all;


  Nr = 800;
  dz_const = 1;
  dz = dz_const*ones(1,Nr);
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);
  Hmax = -zz(end)+dz(end)/2;
  h_shear = 250;

  shear_MITgcm = [0:0.1:2.0]*1e-3;

  for s=1:21
      Shear = shear_MITgcm(s);
      vrelax = h_shear*Shear *(1+ tanh( (zz    +Hmax/2) / (h_shear/2) )) /2;
      mean_shear(s) = (vrelax(450)-vrelax(350)) / (zz(450)-zz(350));
  end

load('MITgcm_growth_linearShear.mat')

figure(1)
clf
set(gcf,'Color','w')
plot(shear_MITgcm,growth_MITgcm,'LineWidth',2)
hold on;
set(gca,'Fontsize',22)
grid on;grid minor;
xlabel('Shear (1/s)')
ylabel('(1/hour)')
title('Growth rate (1/hour)')

load('MITgcm_growth_tanhShear.mat')
hold on;
plot(shear_MITgcm,growth_MITgcm,'LineWidth',2)
plot(mean_shear,growth_MITgcm,'-.','LineWidth',3)

legend('Linear shear','maximum tanh shear','mean tanh shear')


