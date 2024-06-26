
addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo0_H500_s0.0016dz1dx3ln200n-20sm100_kv2e-4';
expdir = '../exps_hires/';
loadexp;

rhoConst = 999.8;

%%% Frequency of diagnostic output
dumpFreq = abs(diag_frequency(1)); 
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);
nDumps = length(dumpIters);

