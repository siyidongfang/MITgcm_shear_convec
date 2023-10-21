%%%
%%% createPBSfile.m
%%%
%%% Writes out a script that runs the MITgcm. Takes as arguments the name
%%% of the experiment, the number of nodes that will be used, and the wall
%%% time in hours (after which the computation will be cut off regardless 
%%% of whether it is complete).
%%%
function createPBSfile_Cheyenne (dirname,expname,nodes,walltime,cluster_path)

  %%% Limit walltime
  % walltime = min(ceil(walltime),12);
  walltime = 12;
  
  %%% Set variables correctly for PBS
  headertext = [...
  '#!/bin/bash \n' ...
  '#PBS -N ',expname,' \n' ...
  '#PBS -A UMIT0054 \n' ...
  '#PBS -m abe \n' ...
  '#PBS -M y_si@mit.edu \n' ...
  '#PBS -q regular \n' ...
  '#PBS -l walltime=',num2str(walltime),':00:00 \n' ...
  '#PBS -l select=',num2str(ceil(nodes/36)),':ncpus=36:mpiprocs=36:ompthreads=1 \n' ...
  '#PBS -o output.txt \n' ...
  '#PBS -e errors.txt \n' ...
  '\n' ...
  'source /etc/profile.d/modules.csh \n' ...
  'module load impi \n' ...
  '\n' ...
  ...%'cd ',cluster_path,'\n' ...
  'mpirun -n ',num2str(nodes),' ./mitgcmuv'];
  
  %%% Open template script
  templatename = './DEFAULTS/results/run_mitgcm';
  tfid = fopen(templatename,'r');
  if (tfid == -1)
    error(['Could not read template PBS run file: ',templatename]);
  end
  
  %%% Open output script and write header text
  wfid = fopen(fullfile(dirname,'run_mitgcm'),'w');
  if (wfid == -1)
    error('Could not open PBS script file');
  end
  fprintf(wfid,headertext); 
  
  %%% Copy over all template text
  count = 0;
  while true
    count = count+1;
    nextline = fgetl(tfid);    
    if ~ischar(nextline)
      break;
    end
    if (count <= 5) %%% skip placeholder comments at start of template
      continue;
    end
    fprintf(wfid,'%s\n',nextline);    
  end
  
  %%% Close files when we're done
  fclose(wfid);
  fclose(tfid);

end

