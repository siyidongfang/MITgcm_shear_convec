
close all;clear;
fontsize = 20;

% expdir = '/Volumes/MIT/MITgcm_shear_convec/instability/experiments/lambda';
% expdir = '/nobackup1/y_si/MITgcm_shear_convec/instability/experiments/lambda';
% expdir = 'exps_test/lores_nu1e-5_lambda';
% lambda_parm = [400 450 550 650 700 750 800 850 1000 1200:200:2400 2800:200:5000 6000 8000:1000:12000];

expdir = 'exps_RK4_nu5e-6/lambda'
% lambda_parm = round(10.^[2:0.1:3.4 3.6 3.8 4])
% lambda_parm = round(10.^[2 2.1 2.4:0.1:3.4 3.6 3.8 4])
lambda_parm = 100
Shear_parm=[0.1:0.2:2.1]*1e-3;


for Nexp_lambda = 1:length(lambda_parm)
    lambda = lambda_parm(Nexp_lambda)

    for Nexp_shear =1:length(Shear_parm)
        Shear = Shear_parm(Nexp_shear);

        expname = ['H300_topo4_Pt43200_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        exppath = [expdir num2str(lambda) '/' expname];
        clear uuu www psi U0 NTtide tt Nr Nt Utide ttd t1hour zz fit_span zzd

        load([exppath '/output.mat'],...
            'www','re_psi','U0','NTtide','tt','Nr','Nt',...
            'ttd','t1hour','zz','zzd','dz')
            uuu = U0*(re_psi(:,2:Nr+1)-re_psi(:,1:Nr))/dz;
            
            fit_span = round(Nt*2/3):Nt-1;

            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
            TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            TPE = 0;
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = mean(KE_PE,2)';
            xxplot = ttd/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(Nexp_lambda,Nexp_shear) = p(1);
            
            
            [y_fit,delta_fit] = polyval(p,xxplot,S);
            h=figure(1);
            clf;set(gcf,'color','w');
            plot(xxplot,yyplot,'LineWidth',2)
            hold on
            plot(xxplot(fit_span),y_fit(fit_span),'LineWidth',2)
            grid on;grid minor;
            set(gca,'Fontsize',fontsize);
            ylim([p(2)-3 p(2)+p(1)*max(xxplot)+2])
            xlabel('$t$','Interpreter','Latex')
            ylabel('$\ln(e)/2$','Interpreter','Latex')
            saveas(h,[exppath 'KE.png'])

       
    end
    
end

save('GrowthRate_RK4_nu5e-6.mat','lambda_parm','Shear_parm','GrowthRate','fit_span')

            
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


