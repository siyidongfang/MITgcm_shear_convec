
% Bethan: For the temperature (insitu temp) and salinity (practical salinity), 
% I’ve attached 7 CTD stations that were done as part of a 25 hour 
% time series near the axis of the canyon about half way along the 
% northern branch, plus 8 CTD stations that were during a transect 
% along the full length of the canyon. 

addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/;
addpath CTD

addpath ../colormaps/
load_colors
load CTD_stations.mat

fontsize = 16;
P_all = NaN.*zeros(1306,15);
pt_all = NaN.*zeros(1306,15);
psal_all = NaN.*zeros(1306,15);
N2_all = NaN.*zeros(1306-1,15);
Pmid_all = NaN.*zeros(1306-1,15);
SA_all = NaN.*zeros(1306,15);
CT_all = NaN.*zeros(1306,15);
lat_mean = 0;

for nf = 1:15
% for nf = 1
    Nz = length(CTD.Stn(nf).pressure);
    pressure = CTD.Stn(nf).pressure;
    insitu_temp = CTD.Stn(nf).insitu_temp;
    psal = CTD.Stn(nf).psal;
    lat = CTD.Stn(nf).lat;
    lon = CTD.Stn(nf).lon;

    % Calculate buoyancy frequency
    [SA, in_ocean] = gsw_SA_from_SP(psal,pressure,lon,lat);
    CT = gsw_CT_from_t(SA,insitu_temp,pressure);
    [N2, p_mid] = gsw_Nsquared(SA,CT,pressure,lat);

    p_ref = 0;
    pt = gsw_pt_from_t(SA,insitu_temp,pressure,p_ref);

    P_all(1:Nz,nf) = pressure;
    pt_all(1:Nz,nf) = pt;
    psal_all(1:Nz,nf) = psal;
    SA_all(1:Nz,nf) = SA;
    CT_all(1:Nz,nf) = CT;
    N2_all(1:Nz-1,nf) = N2;
    Pmid_all(1:Nz-1,nf) = p_mid;

    lat_all(nf) = lat;
    lon_all(nf) = lon;

end

lat_mean7 = mean(lat_all(1:7));

pt_mean7 = mean(pt_all(:,1:7),2,'omitnan');
CT_mean7 = mean(CT_all(:,1:7),2,'omitnan');
psal_mean7 = mean(psal_all(:,1:7),2,'omitnan');
SA_mean7 = mean(SA_all(:,1:7),2,'omitnan');
pp7 = mean(P_all(:,1:7),2,'omitnan');
[N2_mean7, p_mid7] = gsw_Nsquared(SA_mean7,CT_mean7,pp7,lat_mean7);
N2_mean7(N2_mean7<=0)=NaN;
% pt_mean15 = mean(pt_all,2,'omitnan');
% psal_mean15 = mean(psal_all,2,'omitnan');
% SA_mean15 = mean(SA_all,2,'omitnan');
% pp15 = mean(P_all,2,'omitnan');
% [N2_mean15, p_mid15] = gsw_Nsquared(SA_mean15,pt_mean15,pp15,lat_mean);

lat15 = mean([lat_all(8:15) lat_mean7]);
pt_mean15 = mean([pt_all(:,8:15) pt_mean7],2,'omitnan');
psal_mean15 = mean([psal_all(:,8:15) psal_mean7],2,'omitnan');
ref_pres_surf = 0;
SA_mean15 = gsw_SA_from_SP(psal_mean15,ref_pres_surf,-11.9,lat15);  
CT_mean15 = gsw_CT_from_pt(SA_mean15,pt_mean15); 
pp15 = mean(P_all,2,'omitnan');
[N2_mean15, p_mid15] = gsw_Nsquared(SA_mean15,CT_mean15,pp15,lat15);


p1 = p_mid15(find(N2_mean15<0))
% Hmax = 2.4092e+03

negN2idx = find(N2_mean15<0);
negN2idx_model = negN2idx(1:5);

N2_mean15(N2_mean15<=0)=NaN;


for i=1:5
    zidx = negN2idx_model(i)-20:negN2idx_model(i)+20;
    psal_mean15(zidx) = smooth(smooth(smooth(smooth(smooth(smooth(psal_mean15(zidx)))))));
    pt_mean15(zidx) = smooth(smooth(smooth(smooth(smooth(smooth(pt_mean15(zidx)))))));
end
SA_mean15 = gsw_SA_from_SP(psal_mean15,ref_pres_surf,-11.9,lat15);  
CT_mean15 = gsw_CT_from_pt(SA_mean15,pt_mean15); 
[N2_mean15, p_mid15] = gsw_Nsquared(SA_mean15,CT_mean15,pp15,lat15);
p2 = p_mid15(find(N2_mean15<0));


figure(20)
subplot(1,2,1)
plot(pt_mean15,pp15);axis ij;
subplot(1,2,2)
plot(psal_mean15,pp15);axis ij;

save('CTD/CTD.mat','P_all','pt_all','psal_all','SA_all','CT_all','N2_all','Pmid_all','lat_all','lon_all',...
    'lat_mean7','pt_mean7','CT_mean7','psal_mean7','SA_mean7','pp7','N2_mean7','p_mid7',...
    'pt_mean15','psal_mean15','SA_mean15','pp15','N2_mean15','p_mid15','lat15')

%%
% color1=red;color2=orange;color3=pink;color4=yellow;color5=brown;color6=lightblue;color7=purple;color8=green;

figure(1);
clf;set(gcf,'Color','w','Position',[237 206 1117 376])
subplot(1,3,1)
h1 = plot(pt_all(:,1:7),P_all(:,1:7),'--','Color',gray);
axis ij;grid on;grid minor;
hold on;
h2 = plot(pt_all(:,8:15),P_all(:,8:15),'LineWidth',2);
% h8 = plot(pt_all(:,8),P_all(:,8),'Color',color1,'LineWidth',2);
% h9 = plot(pt_all(:,9),P_all(:,9),'Color',color2,'LineWidth',2);
% h10 = plot(pt_all(:,10),P_all(:,10),'Color',color3,'LineWidth',2);
% h11 = plot(pt_all(:,11),P_all(:,11),'Color',color4,'LineWidth',2);
% h12 = plot(pt_all(:,12),P_all(:,12),'Color',color5,'LineWidth',2);
% h13 = plot(pt_all(:,13),P_all(:,13),'Color',color6,'LineWidth',2);
% h14 = plot(pt_all(:,14),P_all(:,14),'Color',color7,'LineWidth',2);
% h15 = plot(pt_all(:,15),P_all(:,15),'Color',color8,'LineWidth',2);
h1_mean = plot(pt_mean15,pp15,'LineWidth',4,'Color','k');hold off;
title('Potential temperature')
xlabel('(degC)');ylabel('Pressure')
set(gca,'Fontsize',fontsize)
ylim([0 2700])
xlim([2 15])
% legend([h1_mean,h2'],'mean','1','2','3','4','5','6','7','8')
subplot(1,3,2)
plot(psal_all(:,1:7),P_all(:,1:7),'--','Color',gray);
axis ij;grid on;grid minor;
hold on;
plot(psal_all(:,8:15),P_all(:,8:15),'LineWidth',2);
plot(psal_mean15,pp15,'LineWidth',4,'Color','k');hold off;
title('Salinity')
ylim([0 2700])
xlabel('(psu)');ylabel('Pressure')
set(gca,'Fontsize',fontsize)
subplot(1,3,3)
plot(log10(N2_all(:,1:7)),Pmid_all(:,1:7),'--','Color',gray);
axis ij;grid on;grid minor;
hold on;plot(log10(N2_all(:,8:15)),Pmid_all(:,8:15))
plot(log10(N2_mean15),p_mid,'LineWidth',2,'Color','k');hold off;
title('log10(N^2)')
ylim([0 2700])
xlim([-9 -3])
xlabel('(s^{-2})');ylabel('Pressure')
set(gca,'Fontsize',fontsize)

% figdir = '/Users/csi/MITgcm_BLT/analysis/NCAR_proposal/';
% print('-djpeg','-r200',[figdir 'CTD.jpeg']);

