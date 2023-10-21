%%%
%%% load_all.m
%%%
%%% A convenient script to make plots for all simulations



    %%% Add path
    addpath /Users/ysi/MITgcm_shear_convec/analysis/;
    addpath /Users/ysi/MITgcm_shear_convec/analysis/functions;
    addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps;
    addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/cmocean/;
    addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/customcolormap/;

    list_exps;
    LinearEOS = true;
    expname = EXPNAME{ne}
    loadexp;
    rhoConst = 999.8;
    load_colors;

    fontsize = 18;
    figdir = [exppath '/img/']; 
    prodir = '/Users/ysi/MITgcm_shear_convec/products/';
    
    dumpFreq = abs(diag_frequency(1)); 
    nDumps = floor(nTimeSteps*deltaT/dumpFreq);
    dumpIters = round((1:nDumps)*dumpFreq/deltaT);
    dumpIters = dumpIters(dumpIters > nIter0);


    % plot_uvtn;
    % plot_RuanFig2;



