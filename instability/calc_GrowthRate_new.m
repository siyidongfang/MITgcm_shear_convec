
close all;
clear;
fontsize = 22;

dirname = 'randn_topo4_0Center';
constants_tanh
expfolder = [dirname '/lambda'];

GrowthRate = NaN.*zeros(length(lambda_parm),length(Shear_parm));
grow =  NaN.*zeros(1,length(Shear_parm));

parfor Nexp_lambda = 1:length(lambda_parm)
    Nexp_lambda
    lambda = lambda_parm(Nexp_lambda);
    grow = zeros(1,length(Shear_parm));

    for Nexp_shear =1:length(Shear_parm)
        % Nexp_shear
        Shear = Shear_parm(Nexp_shear);

        expname = ['topo4_H500_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '_'];
        expdir = [expfolder num2str(lambda) '/' expname];
        % clear uuu www psi NTtide tt Nr Nt Utide tt t1hour zz fit_span

        fname = [expdir 'output.mat'];
        if(isfile(fname))
            
            pb2 =load_func(fname);
            grow(Nexp_shear) = pb2.pb2(1);
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

save(['products/grow_' dirname '.mat'],'a','b','lambda_Floquet','growth_Floquet','shear_Floquet','GrowthRate_Floquet')




function pb2 =load_func(fname)
        pb2 = load( fname ,'pb2');
end
