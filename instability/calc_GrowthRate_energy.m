
close all;clear;
fontsize = 20;

% expdir = '/Volumes/MIT/MITgcm_shear_convec/instability/experiments/lambda';
% expdir = '/nobackup1/y_si/MITgcm_shear_convec/instability/experiments/lambda';
expdir = 'exps_test/lores_lambda';

lambda_parm = [400 450 550 650 700 750 800 850 1000 1200:200:2400 2800:200:5000 6000 8000:1000:12000];
Shear_parm = [0.1:0.1:2.5]*1e-3;


% for Nexp_lambda = 1:length(lambda_parm)
for Nexp_lambda = 1
    lambda = lambda_parm(Nexp_lambda)

    for Nexp_shear =1:length(Shear_parm)
        Shear = Shear_parm(Nexp_shear);

        expname = ['H300_topo4_Pt43200_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        
        clear re_buoy uuu www re_buoyd U0 NTtide tt Nr Nt Utide ttd t1hour zz fit_span zzd

        load([expdir num2str(lambda) '/' expname '/output.mat'],...
            're_buoy','uuu','www','re_buoyd','U0','NTtide','tt','Nr','Nt',...
            'Utide','ttd','t1hour','zz','zzd')
            
            fit_span = round(Nt*0.4):Nt-1;
            
            %%% Option 1
            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
            TKE = 0.25* ( (uuu/U0).^2+(0.5*(www(:,1:Nr)+www(:,2:Nr+1))/U0).^2 );
            TPE =  0.25* ( re_buoy.^2 );
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = mean(KE_PE,2)';
            
            xxplot = tt;
            yyplot = log(KE_PE_zavg)/2;
            [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(1,Nexp_lambda,Nexp_shear) = p(1);
            
            %%% Option 2
            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
            TKE1 = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            TKE2 = uuu.*(Utide*U0);
            TKE = TKE1+TKE2;
            
            TPE = -repmat(zzd,[length(re_buoyd) 1]).*re_buoyd;
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = abs(mean(KE_PE,2))';
            xxplot = ttd/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(2,Nexp_lambda,Nexp_shear) = p(1);
            
            %%% Option 3
            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
            TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            TPE = -repmat(zzd,[length(re_buoyd) 1]).*re_buoyd;
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = abs(mean(KE_PE,2))';
            xxplot = ttd/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(3,Nexp_lambda,Nexp_shear) = p(1);
            
            %%% Option 4
            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
            TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            TPE = abs(-repmat(zzd,[length(re_buoyd) 1]).*re_buoyd);
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = mean(KE_PE,2)';
            xxplot = ttd/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(4,Nexp_lambda,Nexp_shear) = p(1);
            
            
            %%% Option 5
            clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
            TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            TPE = 0;
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = mean(KE_PE,2)';
            xxplot = ttd/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            GrowthRate(5,Nexp_lambda,Nexp_shear) = p(1);
            
            
            [y_fit,delta_fit] = polyval(p,xxplot,S);
            figure(1)
            clf;set(gcf,'color','w');
            plot(xxplot,yyplot,'LineWidth',2)
            hold on
            plot(xxplot(fit_span),y_fit(fit_span),'LineWidth',2)
            grid on;grid minor;
            set(gca,'Fontsize',fontsize);
            ylim([p(2)-3 p(2)+p(1)*max(xxplot)+2])
            xlabel('$t$','Interpreter','Latex')
            ylabel('$\ln(e)/2$','Interpreter','Latex')


        
    end
    
end

save('GrowthRate_lores_H300_Nu2e-4.mat','lambda_parm','Shear_parm','GrowthRate','fit_span')

% save('GrowthRate_lores_H300_largeNuKappa.mat','lambda_parm','Shear_parm','GrowthRate','fit_span')

