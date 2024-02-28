




fit_span = round(Nt*0.4):Nt;

%%% Option 1
clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S
TKE = 0.25* ( (uuu/U0).^2+(0.5*(www(:,1:Nr)+www(:,2:Nr+1))/U0).^2 );
TPE =  0.25* ( re_buoy.^2 );
KE_PE = TKE+TPE;

KE_PE_zavg = mean(KE_PE,2);

xxplot = tt;
yyplot = log(KE_PE_zavg)/2;
[p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
GrowthRate(1,Nexp_lambda,Nexp_shear) = p(1);

%%% Option 2
clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S
TKE1 = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TKE2 = uuu.*(Utide*U0);
TKE = TKE1+TKE2;

TPE = -repmat(zz,[length(re_buoyd) 1]).*re_buoyd;
KE_PE = TKE+TPE;

KE_PE_zavg = abs(mean(KE_PE,2));
xxplot = ttd/t1hour;
yyplot = log(KE_PE_zavg)/2;
[p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
GrowthRate(2,Nexp_lambda,Nexp_shear) = p(1);

%%% Option 3
clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2
TKE1 = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TKE = TKE1;

TPE = -repmat(zz,[length(re_buoyd) 1]).*re_buoyd;
KE_PE = TKE+TPE;

KE_PE_zavg = abs(mean(KE_PE,2));
xxplot = ttd/t1hour;
yyplot = log(KE_PE_zavg)/2;
[p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
GrowthRate(3,Nexp_lambda,Nexp_shear) = p(1);

%%% Option 4
clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2
TKE1 = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TKE = TKE1;

TPE = abs(-repmat(zz,[length(re_buoyd) 1]).*re_buoyd);
KE_PE = TKE+TPE;

KE_PE_zavg = mean(KE_PE,2);
xxplot = ttd/t1hour;
yyplot = log(KE_PE_zavg)/2;
[p,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
GrowthRate(4,Nexp_lambda,Nexp_shear) = p(1);


% [y_fit,delta_fit] = polyval(p,xxplot,S);
% figure(1)
% clf;set(gcf,'color','w');
% plot(xxplot,yyplot,'LineWidth',2)
% hold on
% plot(xxplot(fit_span),y_fit(fit_span),'LineWidth',2)
% grid on;grid minor;
% set(gca,'Fontsize',fontsize);
% ylim([p(2)-3 p(2)+p(1)*max(xxplot)+2])
% xlabel('$t$','Interpreter','Latex')
% ylabel('$\ln(e)/2$','Interpreter','Latex')









