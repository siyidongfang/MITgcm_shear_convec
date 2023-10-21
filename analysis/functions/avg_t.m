%%%
%%% avg_t.m
%%%
%%% Calculates the time average of the output fields from MITgcm runs.
%%%

%%% Read experiment data
clear diag_fields;
clear diag_timePhase;
clear diag_fileNames;
clear diag_frequency;
loadexp;

Nlayers = length(layers_bounds);

savename = [prodir expname '_tavg_5days.mat'];

%%% Frequency of diagnostic output - should match that specified in
%%% data.diagnostics      
dumpFreq = diag_frequency(1);
nDumps = floor(nTimeSteps*deltaT/dumpFreq);
dumpIters = round((1:nDumps)*dumpFreq/deltaT);
dumpIters = dumpIters(dumpIters > nIter0);

%%% Calculate time average for whichever fields are present
 
for m=1:length(diag_fields)
    m
  if (m == 1) %%% N.B. This won't work if the first field isn't one of those listed below
    flag = '';
  else    
    flag = '-append';
  end
  var_name = diag_fields{m};
  % if  (strcmp(var_name(1:2),'La'))  %%% If this is a diagnostic from the LAYERS package, then use Nlayers in the function readIters_days 
  %     if (diag_frequency(m) > 0)
  %         clear avg var_data
  %       var_data = readIters_hours(exppath,var_name,dumpIters,deltaT,tmin,tmax);   
  %       % var_data = squeeze(var_data);
  %       tempStruct.(var_name) = var_data;   
  %       save(savename,'-struct','tempStruct',var_name,flag);
  %     end     
  % else
  if (diag_frequency(m) > 0)
    clear avg var_data
    var_data = readIters_hours(exppath,var_name,dumpIters,deltaT,tmin,tmax); 
    % var_data = squeeze(var_data);
    tempStruct.(var_name) = var_data;   
    save(savename,'-struct','tempStruct',var_name,flag);
  end
  % end
end

