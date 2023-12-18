
clear;close all;

 figure(1)
 clf;set(gcf,'Color','w','Position',[117 426 872 496])

for ne = 1:4
    load_all
    filename = [expdir expname '/RMSE.mat'];
    load(filename)
   
    figure(1)
    hold on;
    plot(time_h,div_uu_zavg,'LineWidth',2);
end

    set(gca,'Fontsize',fontsize)
    xlabel('Time (hours)')
    % title('Temperature RMSE averaged over the bottom shear layer')
    title('Velocity u RMSE averaged over the bottom shear layer')
    % ylabel('(degC)')
    ylabel('(m/s)')
    grid on;grid minor;
    legend('\Lambda = 1.0\times10^{-3} s^{-1}',...
        '\Lambda = 1.1\times10^{-3} s^{-1}',...
        '\Lambda = 1.2\times10^{-3} s^{-1}',...
        '\Lambda = 1.3\times10^{-3} s^{-1}','Fontsize',fontsize,'Position', [0.1711 0.6063 0.2024 0.2288])