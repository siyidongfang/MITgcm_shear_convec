%%%
%%% readIters_hours.m
%%%
%%% Reads ans sums all iterations of a specified MITgcm output field
%%% between specified times, and calculates the time average.
%%%
function avg = readIters_hours (exppath,field,dumpIters,deltaT,tmin,tmax)
 
  avg = zeros(size(rdmdsWrapper(fullfile(exppath,'results',field),dumpIters(1))));
  navg = 0;
  
  %%% Loop through output iterations
  for n=1:length(dumpIters)
     
    thours =  dumpIters(n)*deltaT/3600;
    
    if ((thours >= tmin) && (thours <= tmax))
     
      avg = avg + rdmdsWrapper(fullfile(exppath,'results',field),dumpIters(n));      
      navg = navg + 1;

    end
    
  end
  % navg
  %%% Calculate average
  if (navg > 0)
    avg = avg / navg;
  else
    error('No output files found');
  end
 
end

