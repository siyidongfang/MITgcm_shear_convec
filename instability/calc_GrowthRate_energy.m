
close all;
clear;
fontsize = 22;

dirname = 'exps_linear_dz0.5';
expfolder = [dirname '/lambda'];
Shear_parm = ([0:0.1:2.0])*1e-3;
lambda_parm = [50 75:25:1000 1050:50:2000 2200:200:3000 3500 4000 4500 5000 7500 10000]; 
% lambda_parm = [round(10.^[1.7:0.05:3 3.1:0.1:3.4 3.6 3.8 4]/10)*10];
lambda_parm = flip(lambda_parm);
lambda_parm = [lambda_parm round(10.^[1.6:-0.1:0.5])];

GrowthRate = NaN.*zeros(length(lambda_parm),length(Shear_parm));
grow =  NaN.*zeros(1,length(Shear_parm));

% for Nexp_lambda = 8:length(lambda_parm)
for Nexp_lambda=47:60
    Nexp_lambda
    lambda = lambda_parm(Nexp_lambda);
    grow = zeros(1,length(Shear_parm));

    for Nexp_shear =1:length(Shear_parm)
        Nexp_shear
        Shear = Shear_parm(Nexp_shear);

        expname = ['topo0_H250_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        expdir = [expfolder num2str(lambda) '/' expname];
        % clear uuu www psi NTtide tt Nr Nt Utide tt t1hour zz fit_span

        fname = [expdir 'output.mat'];
        if(isfile(fname))
            
            load([expdir 'output.mat'])
            make_figures
            % [www re_psi NTtide Nr Nt tt t1hour zz dz re_buoy nu kappa] =load_func(fname);
            uuu = (re_psi(:,2:Nr+1)-re_psi(:,1:Nr))/dz;
            
            fit_span = round(Nt/NTtide*20):Nt-1;

            % clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 pp S 
            TKE = 1/4*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
            Pr = nu/kappa;
            TPE = Pr*re_buoy.^2/4;
            KE_PE = TKE+TPE;
            
            KE_PE_zavg = mean(KE_PE,2,'omitnan');
            xxplot = tt/t1hour;
            yyplot = log(KE_PE_zavg)/2;
            [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
            grow(Nexp_shear) = pp(1);
        end


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

    GrowthRate(Nexp_lambda,:) = grow;

end


%%
GrowthRate (GrowthRate==0)=NaN;
GrowthRate_Floquet = GrowthRate;
growth_Floquet = max(GrowthRate);
shear_Floquet = Shear_parm;
lambda_Floquet = lambda_parm;

figure(3)
set(gcf,'color','w')
plot(Shear_parm,growth_Floquet,'LineWidth',2);
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('Shear (1/s)')
ylabel('Growth rate (1/hour)')

figure(4)
set(gcf,'color','w')
pcolor(shear_Floquet,log10(lambda_Floquet),GrowthRate_Floquet);shading flat;colorbar;
clim([-0.3 0.3]);colormap(redblue)
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('Shear (1/s)')
title('Growth rate (1/hour)')
ylabel('log_{10}(\lambda_x) (m)')

GrowthRate_Floquet(isnan(GrowthRate_Floquet)) = 0;
for ns = 1:21
    [a(ns) b(ns)] = max(GrowthRate_Floquet(:,ns));
end

GrowthRate_Floquet(GrowthRate_Floquet==0)=NaN;

save(['GrowthRate_' dirname '.mat'],'a','b','lambda_Floquet','growth_Floquet','shear_Floquet','GrowthRate_Floquet')




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


function [www re_psi NTtide Nr Nt tt t1hour zz dz re_buoy nu kappa] =load_func(file)
        S = load( file );
        www = S.www;
        re_psi = S.re_psi;
        NTtide = S.NTtide;
        Nr = S.Nr;
        Nt = S.Nt;
        tt = S.tt;
        t1hour = S.t1hour;
        zz = S.zz;
        dz = S.dz;
        re_buoy = S.re_buoy;
        nu = S.nu;
        kappa = S.kappa;
end
