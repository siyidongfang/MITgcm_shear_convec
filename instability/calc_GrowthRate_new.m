
close all;
clear;
fontsize = 22;

dirname = 'new_topo4_linear';
expfolder = [dirname '/lambda'];
% Shear_parm = ([0.1:0.1:2.8])*1e-3;    %%% flat
Shear_parm = ([0.1:0.1:2.0 2.07])*1e-3; %%% topo4
lambda_parm = [50:25:700 750:50:1000 1200:200:2000 2400:400:3200 4000:1000:8000 10000:2000:12000]; 
lambda_parm = flip(lambda_parm);
lambda_parm = [lambda_parm round(10.^[1.6:-0.1:0.7]) 3];

GrowthRate = NaN.*zeros(length(lambda_parm),length(Shear_parm));
grow =  NaN.*zeros(1,length(Shear_parm));

parfor Nexp_lambda = 1:length(lambda_parm)
    Nexp_lambda
    lambda = lambda_parm(Nexp_lambda);
    grow = zeros(1,length(Shear_parm));

    for Nexp_shear =1:length(Shear_parm)
        Nexp_shear
        Shear = Shear_parm(Nexp_shear);

        expname = ['topo4_H250_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        expdir = [expfolder num2str(lambda) '/' expname];
        % clear uuu www psi NTtide tt Nr Nt Utide tt t1hour zz fit_span

        fname = [expdir 'output_new.mat'];
        if(isfile(fname))
            
            pb2 =load_func(fname);
            grow(Nexp_shear) = pb2;
        end



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
for ns = 1:length(shear_Floquet)
    [a(ns) b(ns)] = max(GrowthRate_Floquet(:,ns));
end

GrowthRate_Floquet(GrowthRate_Floquet==0)=NaN;

save(['GrowthRate_' dirname '_new.mat'],'a','b','lambda_Floquet','growth_Floquet','shear_Floquet','GrowthRate_Floquet')




function pb2 =load_func(fname)
        pb2 = load( fname ,'pb2');
end
