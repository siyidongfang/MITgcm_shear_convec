
close all;clear;
fontsize = 20;

expdir = '/Volumes/MIT/MITgcm_shear_convec/instability/experiments/lambda';
lambda_parm = [1000 1200 1400 1600 10000 11000 12000 13000 14000 15000 16000 100000];
Shear_parm = [0.1:0.1:0.8]*1e-3;


for Nexp_lambda = 1:length(lambda_parm)
    lambda = lambda_parm(Nexp_lambda)

    for Nexp_shear =1:length(Shear_parm)
        Shear = Shear_parm(Nexp_shear);

        expname = ['H300_topo4_Pt43200_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        
        load([expdir num2str(lambda) '/' expname '/output.mat'],...
            're_buoy','uuu','www','re_buoyd','U0','NTtide','tt','Nr','Nt',...
            'Utide','ttd','t1hour','zz')

        calc_growth;
        
    end
    
end

save('GrowthRate.mat','lambda_parm','Shear_parm','GrowthRate','fit_span')

