
clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo4_H500_smo100m_s0.0014_dz1dx3ln200n-20'
expdir = '/Volumes/MIT/MITgcm_shear_convec/exps_topo4/';
loadexp;
rhoConst = 999.8;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

botZ =zz(end);
