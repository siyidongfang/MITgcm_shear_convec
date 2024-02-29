
close all;clear;
fontsize = 20;

% expdir = '/Volumes/MIT/MITgcm_shear_convec/instability/experiments/lambda';
expdir = '/nobackup1/y_si/MITgcm_shear_convec/instability/experiments/lambda';
lambda_parm = [400 450 550 650 700 750 800 850 1000 1200:200:2400 2800:200:5000 6000 8000:1000:12000];
Shear_parm = [0.1:0.1:2.0]*1e-3;


for Nexp_lambda = 1:length(lambda_parm)
    lambda = lambda_parm(Nexp_lambda)

    for Nexp_shear =1:length(Shear_parm)
        Shear = Shear_parm(Nexp_shear);

        expname = ['H300_topo4_Pt43200_N0.001_S' num2str(Shear) '_lambda' num2str(lambda) '/'];
        
        clear re_buoy uuu www re_buoyd U0 NTtide tt Nr Nt Utide ttd t1hour zz fit_span zzd

        load([expdir num2str(lambda) '/' expname '/output.mat'],...
            're_buoy','uuu','www','re_buoyd','U0','NTtide','tt','Nr','Nt',...
            'Utide','ttd','t1hour','zz','zzd')

        calc_growth;

        
    end
    
end


save('GrowthRate_new.mat','lambda_parm','Shear_parm','GrowthRate','fit_span')

