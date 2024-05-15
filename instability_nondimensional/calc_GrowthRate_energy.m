
close all;
clear;
fontsize = 22;

% expdir = '/Volumes/MIT/MITgcm_shear_convec/instability/experiments/lambda';
% expdir = '/nobackup1/y_si/MITgcm_shear_convec/instability/experiments/lambda';
% expdir = 'exps_test/lores_nu1e-5_lambda';
% lambda_parm = [400 450 550 650 700 750 800 850 1000 1200:200:2400 2800:200:5000 6000 8000:1000:12000];


expdir = 'exps_test_Nr100/lambda'
lambda_parm = [80 100 130 160 250 320 500 630 1000 1260 2000 2510 6310 10000]
Shear_parm=[0.1:0.2:2.1]*1e-3;


for Nexp_lambda = 1:length(lambda_parm)
    Nexp_lambda
    lambda = lambda_parm(Nexp_lambda);

    for Nexp_shear =1:length(Shear_parm)
        % Nexp_shear
        Shear = Shear_parm(Nexp_shear);

        expname = ['H250_topo4_Pt43200_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        exppath = [expdir num2str(lambda) '/' expname];
        clear uuu www psi U0 NTtide tt Nr Nt Utide ttd t1hour zz fit_span zzd

        load([exppath '/output.mat'],...
            'www','re_psi','U0','NTtide','tt','Nr','Nt',...
            'ttd','t1hour','zz','zzd','dz','re_buoy','nu','kappa')
            uuu = U0*(re_psi(:,2:Nr+1)-re_psi(:,1:Nr))/dz;
            
            fit_span = round(Nt/NTtide*3):Nt-1;

            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 pp S 
            TKE = 1/4*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            TKE = TKE/(0.5*U0^2);
            Pr = nu/kappa;
            TPE = Pr*re_buoy.^2/4;
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = mean(KE_PE,2)';
            xxplot = ttd/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(Nexp_lambda,Nexp_shear) = pp(1);
            % % % pp(1);
            % % % [y_fit,delta_fit] = polyval(pp,xxplot,S);
            % % % 
            % % % % b2 = mean(re_buoy.^2,2)';
            % % % % yyplot_b2 = log(b2)/2;
            % % % % [pb2,S_b2] = polyfit(xxplot(fit_span),yyplot_b2(fit_span),1); 
            % % % % GrowthRate_b2(Nexp_lambda,Nexp_shear) = pb2(1);
            % % % % growth(Nexp_lambda,Nexp_shear) = pb2(1)
            % % % % [y_fit_b2,delta_fit_b2] = polyval(pb2,xxplot,S_b2);
            % % % 
            % % % 
            % % % h=figure(1);
            % % % clf;
            % % % set(gcf,'color','w','Position',[85 222 979 420]);
            % % % plot(xxplot/12,yyplot,'LineWidth',2)
            % % % hold on
            % % % % plot(xxplot/12,yyplot_b2,'LineWidth',2)
            % % % plot(xxplot(fit_span)/12,y_fit(fit_span),':','LineWidth',1.5)
            % % % % plot(xxplot(fit_span)/12,y_fit_b2(fit_span),':','LineWidth',1.5)
            % % % grid on;grid minor;
            % % % set(gca,'Fontsize',fontsize);
            % % % ylim([pp(2)-3 pp(2)+pp(1)*max(xxplot)+2])
            % % % xlabel('$t$ (tidal cycle)','Interpreter','Latex')
            % % % ylabel('$\ln(e)/2$','Interpreter','Latex')
            % % % hold off;
            % % % % legend('TKE','b^2','Position',[0.8141 0.1988 0.0684 0.1393])
            % % % saveas(h,[exppath 'energy.png'])

       
    end
    
end


GrowthRate_Floquet = GrowthRate;
growth_Floquet = max(GrowthRate);
shear_Floquet = Shear_parm;
lambda_Floquet = lambda_parm;

figure(3)
set(gcf,'color','w')
plot(Shear_parm,growth_Floquet,'LineWidth',2);
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('Shear (m/s)')
ylabel('Growth rate (1/hour)')

save('GrowthRate_Floquet_topo4.mat','lambda_Floquet','growth_Floquet','shear_Floquet','GrowthRate_Floquet')

            
% %%% Option 2
% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
% TKE = 0.25* ( (uuu/U0).^2+(0.5*(www(:,1:Nr)+www(:,2:Nr+1))/U0).^2 );
% TPE =  0.25* ( re_buoy.^2 );
% KE_PE = TKE+TPE;
% 
% KE_PE_zavg = mean(KE_PE,2)';
% 
% xxplot = tt;
% yyplot = log(KE_PE_zavg)/2;
% [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% GrowthRate(2,Nexp_lambda,Nexp_shear) = p(1);
% 
% %%% Option 3
% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
% TKE1 = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
% TKE2 = uuu.*(Utide*U0);
% TKE = TKE1+TKE2;
% 
% TPE = -repmat(zzd,[length(re_buoyd) 1]).*re_buoyd;
% KE_PE = TKE+TPE;
% 
% KE_PE_zavg = abs(mean(KE_PE,2))';
% xxplot = ttd/t1hour;
% yyplot = log(KE_PE_zavg)/2;
% [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% GrowthRate(3,Nexp_lambda,Nexp_shear) = p(1);
% 
% %%% Option 4
% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
% TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
% TPE = -repmat(zzd,[length(re_buoyd) 1]).*re_buoyd;
% KE_PE = TKE+TPE;
% 
% KE_PE_zavg = abs(mean(KE_PE,2))';
% xxplot = ttd/t1hour;
% yyplot = log(KE_PE_zavg)/2;
% [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% GrowthRate(4,Nexp_lambda,Nexp_shear) = p(1);
% 
% %%% Option 5
% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
% TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
% TPE = abs(-repmat(zzd,[length(re_buoyd) 1]).*re_buoyd);
% KE_PE = TKE+TPE;
% 
% KE_PE_zavg = mean(KE_PE,2)';
% xxplot = ttd/t1hour;
% yyplot = log(KE_PE_zavg)/2;
% [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% GrowthRate(5,Nexp_lambda,Nexp_shear) = p(1);


