%%%
%%% avg_t_tides.m
%%%
%%% Calculates the hourly average of the output fields from MITgcm runs,
%%% averaged over 10 tidal cycles
%%%

%%% Read experiment data
clear diag_fields;
clear diag_timePhase;
clear diag_fileNames;
clear diag_frequency;
loadexp;

% Nlayers = length(layers_bounds);

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics      
dumpFreq = diag_frequency(1);
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

for h_tides = 1:12

    h_tides
    savename = [prodir expname '_tavg_5days_h' num2str(h_tides) '.mat'];
    %%% Calculate time average for whichever fields are present

    dumpIters_htides = dumpIters(h_tides:12:length(dumpIters));
    % dumpIters_htides = dumpIters(h_tides:12:12);

    for m=1:length(diag_fields)
      if (m == 1) %%% N.B. This won't work if the first field isn't one of those listed below
        flag = '';
      else    
        flag = '-append';
      end
      var_name = diag_fields{m};
      if (diag_frequency(m) > 0)
          clear avg var_data
        var_data = readIters_tides(exppath,var_name,dumpIters_htides);   
        var_data = squeeze(var_data);
        tempStruct.(var_name) = var_data;   
        save(savename,'-struct','tempStruct',var_name,flag);
      end     
    end

end

