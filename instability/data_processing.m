

clear all;
close all;

Shear_parm = ([0:0.1:2.0])*1e-3;
lambda_parm = [round(10.^[1.7:0.05:3 3.1:0.1:3.4 3.6 3.8 4]/10)*10];
lambda_parm = flip(lambda_parm);
lambda_parm = [lambda_parm round(10.^[1.6:-0.1:0.5])];

exppath = 'exps_linear/';

topo = 0;
Hmax = 250;
N = 1e-3;

for Nexp_lambda =40:length(lambda_parm)
    lambda = lambda_parm(Nexp_lambda)
    expfolder = [exppath 'lambda' num2str(lambda) '/'];

    for Nexp_shear =1:length(Shear_parm)
	Nexp_shear
        Shear = Shear_parm(Nexp_shear);
        expdir = [expfolder 'topo' num2str(topo) '_H' num2str(Hmax) ...
            '_N' num2str(N) '_S' num2str(Shear) ...
            '_lambda' num2str(lambda) '/'];
        outputname = [expdir 'output.mat'];
        
        if(isfile(outputname))
            load(outputname)
            clear b0 b_wgrid b0_wgrid b_2 b_3 b_4 p0 p0_ugrid psi psi0 sol1 solinit ...
                z0 z_2 z_3 z_4 zeta dbdz ...
                d2bdz2 d2psidz2 d2zetadz2  dpsidz dUtidedz dzetadt ...
                bq1 bq2 bq3 bq4 bq5 zq1 zq2 zq3 zq4 ...
                k_1b k_1z k_2b k_2z k_3b k_3z k_4b k_4z h ...
                buoy re_zq1 re_zq2 re_zq3 re_zq4 re_bq1 re_bq2 re_bq3 re_bq4 re_bq5 Utide ...
                uuu re_psid re_zetad re_buoyd re_dbdz re_d2bdz2 re_d2zetadz2 
            
            outputname2 = [expdir 'output2.mat'];
            save(outputname2)
        end

    end
    
end



