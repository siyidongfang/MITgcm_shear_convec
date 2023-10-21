%%%
%%% readIters_tides.m
%%%
%%% Reads ans sums all iterations of a specified MITgcm output field
%%% between specified times, and calculates the time average.
%%%
function avg = readIters_tides (exppath,field,dumpIters_htides)
 
  avg = zeros(size(rdmdsWrapper(fullfile(exppath,'results',field),dumpIters_htides(1))));
  navg = 0;
  
  %%% Loop through output iterations
  for n=1:length(dumpIters_htides)

      avg = avg + rdmdsWrapper(fullfile(exppath,'results',field),dumpIters_htides(n));      
      navg = navg + 1;

  end
  
  %%% Calculate average
  if (navg > 0)
    avg = avg / navg;
  else
    error('No output files found');
  end
 
end

