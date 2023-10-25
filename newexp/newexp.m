
%%%
%%% newexp.m
%%%
%%% Convenience script to set up folders and options for a new MITgcm run
%%% in a way that allows:
%%% - automatic generation of basic folder and file structure
%%% - automatic calculation of inter-dependent input parameters
%%% - compatibility between required numbers of processors
%%%
%%% Sets up a new experiment folder with subdirectories for the build, the 
%%% code opti ons, the input files, and a results folder. Creates a 'build'
%%% script in the build folder and a 'run' script in the results folder.
%%% Generates the SIZE.h file in the code folder based on parameters
%%% specified here, and copies other code files from the DEFAULTS/code 
%%% folder. Generates all 'eedata' and some 'data' parameters in the input 
%%% folder based on parameters specified here and code in create_data 
%%% function. Other parameters are copied from the DEFAULTS/input folder.
%%%
%%% NOTE: 'expname' MUST NOT be set to 'DEFAULTS'
%%%

function newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)
 
  select_DEFAULTS = 'DEFAULTS_new';
   
%    addpath /data/MITgcm_ASF-csi/newexp_utils/;
  addpath ../newexp_utils/;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% USER-SPECIFIED OPTIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%clo%%%%%%%%%%  
  
  
  %%% Local directory in which to create experiments
  foldername = 'MITgcm_shear_convec';
  runsdir = ['/Users/ysi/' foldername '/'];

  %%% Experiment subdirectories 
  builddir = 'build';
  codedir = 'code';
  inputdir = 'input';
  resultsdir = 'results';
  imgdir = 'img';
 
  %%% List terminator character for parameter files - may be '/' or '&'
  %%% depending on operating system
  listterm = '&';

  %%% Line feed character - important for .sh shell script files
  %%% On unix this probably need to be '\n', in windows '\r\n\'
  lf = '\n';    
  
  nSx = 1; %%% no. of tiles per processor in x-direction
  nSy = 1; %%% no. of tiles per processor in y-direction
  nTx = 1; %%% no. of threads per processor in x-direction
  nTy = 1; %%% no. of threads per processor in y-direction
  OLx = 4; %%% no. of overlapping x-gridpoints per tile
  OLy = 4; %%% no. of overlapping y-gridpoints per tile    

  % %%% Set-up for Derecho, 64 CPUs per nodes
  % opt_file = 'linux_amd64_ifort+impi_derecho';
  % use_mpi = true; %%% set true for parallel processing
  % use_pbs = true; %%% set true for execution via PBS
  % cluster = 'derecho';    
  % nPx = 1; %%% no. of processors in x-direction
  % nPy = 250; %%% no. of processors in y-direction
  % sNx = 1; %%% no. of x-gridpoints per tile
  % sNy = Ny/nPy; %%% no. of y-gridpoints per tile
 
  % %%% Set-up for Cheyenne, 36 CPUs per nodes
  % opt_file = 'linux_amd64_ifort+impi'; %%% options file name
  % use_mpi = true; %%% set true for parallel processing
  % use_pbs = true; %%% set true for execution via PBS
  % cluster = 'cheyenne';
  % nPx = 1; %%% no. of processors in x-direction
  % nPy = 216; %%% no. of processors in y-direction
  % sNx = 1; %%% no. of x-gridpoints per tile
  % sNy = Ny/nPy; %%% no. of y-gridpoints per tile

  %%% Set-up for Engaging
  opt_file = 'linux_amd64_ifort+impi_engaging'; %%% options file name
  use_mpi = true; %%% set true for parallel processing
  use_pbs = true; %%% set true for execution via PBS
  cluster = 'engaging';    
  nPx = 100; %%% no. of processors in x-direction
  nPy = 1; %%% no. of processors in y-direction
  sNx = Nx/nPx; %%% no. of x-gridpoints per tile
  sNy = 1; %%% no. of x-gridpoints per tile
  % sNy = Ny/nPy; %%% no. of y-gridpoints per tile


  % %%% Set-up for Si's Mac Pro - barotropic test case
  % opt_file = 'darwin_arm64_gfortran'; %%% options file name
  % use_mpi = false; %%% set true for parallel processing
  % use_pbs = false; %%% set true for execution via PBS
  % cluster = 'none';    
  % sNx = 1; %%% no. of x-gridpoints per tile
  % sNy = 75; %%% no. of y-gridpoints per tile
  % nPx = 1; %%% no. of processors in x-direction
  % nPy = 1; %%% no. of processors in y-direction
  % Nr = 66; %%% no. of z-gridpoints

  % %%% Set-up for Hoffman2
  % opt_file = 'hoffman2_ucla'; %%% options file name
  % use_mpi = true; %%% set true for parallel processing
  % use_pbs = true; %%% set true for execution via PBS
  % cluster = 'hoffman2';
  % queue = 'all.q';
  % if (is_hires) 
  %   sNx = 1; %%% no. of x-gridpoints per tile
  %   sNy = 10; %%% no. of y-gridpoints per tile
  %   nPx = 1; %%% no. of processors in x-direction
  %   nPy = 50; %%% no. of processors in y-direction
  %   Nr = 50; %%% no. of z-gridpoints
  % else
  %   sNx = 1; %%% no. of x-gridpoints per tile
  %   sNy = 8; %%% no. of y-gridpoints per tile
  %   nPx = 1; %%% no. of processors in x-direction
  %   nPy = 25; %%% no. of processors in y-direction
  %   Nr = 50; %%% no. of z-gridpoints
  % end


  %%% Uploading/downloading parameters 
  switch(cluster)
    case 'derecho'
        
      username = 'ysi';
      clustername = 'derecho.hpc.ucar.edu';
      toolsdir = fullfile(['/glade/derecho/scratch/',username,'/' foldername '/tools/']);
      clusterdir = fullfile(['/glade/derecho/scratch/',username,'/' foldername '/'],batch_name);      

    case 'cheyenne'
        
      username = 'ysi';
      clustername = 'cheyenne.ucar.edu';
      toolsdir = fullfile(['/glade/scratch/',username,'/' foldername '/tools/']);
      clusterdir = fullfile(['/glade/scratch/',username,'/' foldername '/'],batch_name);      
          
    case 'engaging'
        
      username = 'y_si';
      clustername = 'eofe7.mit.edu';
      toolsdir = fullfile(['/nobackup1c/users/',username,'/' foldername '/tools/']);
      clusterdir = fullfile(['/nobackup1c/users/',username,'/' foldername '/'],batch_name);      
        
    case 'gordon'
        
      username = 'ysi';
      clustername = 'gordon.sdsc.xsede.org';
      toolsdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/tools/');
      clusterdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/',batch_name);      
           
    case 'stampede'
        
      username = 'tg854737';
      clustername = 'stampede2.tacc.utexas.edu';
      toolsdir = fullfile('/scratch/06174/tg854737/MITgcm_BLT_ysi/tools/');
      clusterdir = fullfile('/scratch/06174/tg854737/MITgcm_BLT_ysi/',batch_name);      
      
    case 'comet'
        
      username = 'csi';
      clustername = 'comet.sdsc.xsede.org';
      toolsdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/tools/');
      clusterdir = fullfile('/oasis/projects/nsf/',acct,username,'/MITgcm_WS/experiments/',batch_name);      
            
    case 'ardbeg' 
  
      username = 'csi';
      clustername = 'caolila.atmos.ucla.edu';
      toolsdir = '/data2/csi/MITgcm_BLT_ysi/tools/';
      clusterdir = fullfile('/data2/csi/MITgcm_BLT_ysi/',batch_name);      
      
    case 'hoffman2'
        
      username = 'csi';
      clustername = 'hoffman2.idre.ucla.edu';
      toolsdir = '/u/scratch/c/csi/MITgcm_BLT_ysi/tools/';
      clusterdir = ['/u/scratch/c/csi/MITgcm_BLT_ysi/',batch_name];
      
    otherwise %%% Defaults to 'none'
       
      username = 'ysi';
      clustername = '';
      toolsdir = '/Users/ysi/MITgcm_BLT/tools/';
      clusterdir = fullfile('/Users/ysi/MITgcm_BLT/exps_test/',batch_name);  
     
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NON-USER-SPECIFIED PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Paths to sub-directories
  dirpath = fullfile(runsdir,batch_name);
  exppath = fullfile(dirpath,exp_name);
  buildpath = fullfile(exppath,builddir);
  codepath =  fullfile(exppath,codedir);
  inputpath = fullfile(exppath,inputdir);
  resultspath = fullfile(exppath,resultsdir);
  imgpath = fullfile(exppath,imgdir);
  %%% We have to use MPI if we're using PBS
  if (use_pbs)
    use_mpi = true;    
  end

  %%% If use_mpi is false then we can only have one processor
  if ((use_mpi == false) && ((nPx~=1) || (nPy~=1)))
    error('Only one processor allowed for non-MPI computation');
  end
  
  %%% Calculate total grid size and number of nodes
  Nx = sNx*nSx*nPx;
  Ny = sNy*nSy*nPy;
  nodes = nPx*nPy;



  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIRECTORIES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%


  %%% Open experiment top-level directory
  [dir_success,dir_msg,dir_msgid] = mkdir(exppath);
  if (dir_success == 0)
    error(strcat(['Could not open ',exp_name,' : ',num2str(dir_msgid),' : ',dir_msg]));
  end

  %%% Open sub-directories
  subdirnames = {builddir,codedir,inputdir,resultsdir,imgdir};
  for n=1:1:length(subdirnames)     
    [dir_success,dir_msg,dir_msgid] = mkdir(exppath,subdirnames{n});
    if (dir_success == 0)
      error(strcat(['Could not open ',exppath,subdirnames{n},' : ',num2str(dir_msgid),' : ',dir_msg]));
    end
  end


  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%

  
  %%% Copy other files across
  codelist = dir(['./' select_DEFAULTS '/input/']);
  for n=1:1:length(codelist)
    %%% Ignore hidden files
    if (codelist(n).name(1) == '.')
      continue;
    end    
    copyfile(fullfile(['./' select_DEFAULTS '/input/'],codelist(n).name),fullfile(inputpath,codelist(n).name));
  end   


  [nTimeSteps,h,tNorth,sNorth,rho_north,N]...
    = setParams(exp_name,inputpath,codepath,imgpath,listterm,Nx,Ny,Nr,Atide,randtopog_height,randtopog_length,run_type);


  %%% Generate 'eedata'
  create_eedata(inputpath,listterm,nTx,nTy);     

  
  %%%%%%%%%%%%%%%%%%%%%%
  %%%%% CODE FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%


  %%% Generate SIZE.h and just copy other code files
  createSIZEh(codepath,sNx,sNy,nSx,nSy,nPx,nPy,OLx,OLy,Nr);
  codelist = dir(['./' select_DEFAULTS '/code/']);
  for n=1:1:length(codelist)
    if (codelist(n).name(1) == '.')
      continue;
    end
    copyfile(fullfile(['./' select_DEFAULTS '/code/'],codelist(n).name),fullfile(codepath,codelist(n).name));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ESTIMATE WALL TIME %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Computation time (in hours) per gridpoint (in space and time) 
  %%% assigned to each processor.
  %%% Estimated for a single Fram core.
  switch (cluster)
    case 'derecho'
      alpha = 1.59e-9*1.844;
    case 'cheyenne'
      alpha = 1.59e-9*6;
    case 'engaging'
      alpha = 1.59e-9*3;
    case 'gordon'
      alpha = 0.63e-9;
    case 'stampede'
      alpha = 1.59e-9;
    case 'comet'
      alpha = 0.63e-9;
    case 'hoffman2'
      alpha = 0.63e-9;
    otherwise %%% Defaults to Ardbeg
      alpha = 0.63e-9;
  end  
  
  %%% Estimated total computation time in hours (assumes one tile per
  %%% processor). 
  
  %%% This tends to overestimate the computation time when OLx is
  %%% comparable to sNx or OLy is comparable to sNy.
  % comptime = alpha*(sNx+2*OLx)*(sNy+2*OLy)*Nr*Nt  
  
  %%% This seems to provide a decent estimate when OLx is
  %%% comparable to sNx or OLy is comparable to sNy; 'ghost' gridpoints
  %%% require perhaps half as much processing as 'real' gridpoints.
  comptime = alpha*(sNx+OLx)*(sNy+OLy)*Nr*nTimeSteps
    
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR BUILDING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Build commands - depend on whether MPI is used
  if (use_mpi)
    mpistr = '-mpi ';
  else 
    mpistr = '';
  end
  if(~strcmp(cluster,'engaging'))
  buildcommands = strcat([... 
    'rm my_opt_file',lf,...
    'ln -s ',fullfile(toolsdir,'build_options',opt_file),' my_opt_file',lf, ...
    fullfile(toolsdir,'genmake2'),' ',mpistr,'-mods=../code -of=my_opt_file',lf, ...
    'make depend',lf, ...
    'make --always-make -j 2',lf,]);
  else
  buildcommands = strcat([... 
    'source /home/y_si/.bashrc',lf,...
    'rm my_opt_file',lf,...
    'ln -s ',fullfile(toolsdir,'build_options',opt_file),' my_opt_file',lf, ...
    fullfile(toolsdir,'genmake2'),' ',mpistr,'-mods=../code -of=my_opt_file',lf, ...
    'make depend',lf, ...
    'make --always-make -j 2',lf,]);
  end

  %%% Create the 'build' shell file
  fid = fopen(fullfile(buildpath,'build.sh'),'w');
  if (fid == -1)
    error('Could not open build script for writing');
  end
  fprintf(fid,buildcommands);
  fclose(fid);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR RUNNING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cl%%%%%%%%%%%


  %%% Commands to link input and build folders to results folder
  runcommands = [...
    'ln -s ../',inputdir,'/* ./ ',lf, ...
    'ln -s ../',builddir,'/mitgcmuv ',lf];

  %%% Execution command depends on whether MPI is used
  if (use_pbs)
    switch (cluster)
       case 'derecho'
        createPBSfile_Derecho(resultspath,exp_name,nodes,2*comptime,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
      case 'cheyenne'
        createPBSfile_Cheyenne(resultspath,exp_name,nodes,2*comptime,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
      case 'engaging'
        createPBSfile_Engaging(resultspath,exp_name,nodes,2*comptime,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
      case 'gordon'
        createPBSfile_Gordon(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
      case 'stampede'
        createPBSfile_Stampede(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));       
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
      case 'stampede2'
        createPBSfile_Stampede2(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));       
        runcommands = [runcommands,'sbatch run_mitgcm',lf];
      case 'comet'
        createPBSfile_Comet(resultspath,exp_name,nodes,2*comptime,acct,fullfile(clusterdir,exp_name,resultsdir));        
        runcommands = [runcommands,'sbatch run_mitgcm_comet',lf];
      case 'hoffman2'
        createPBSfile_Hoffman(resultspath,exp_name,nodes);        
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
      otherwise %%% Defaults to Ardbeg
        createPBSfile(resultspath,exp_name,nodes,queue);
        runcommands = [runcommands,'qsub run_mitgcm > output.txt',lf];
    end    
  else
    if (use_mpi)
      runcommands = [runcommands,'mpirun -np ',num2str(nodes), ...
                      ' ./mitgcmuv > output.txt',lf];
    else
      runcommands = [runcommands,'./mitgcmuv > output.txt',lf];
    end
  end

  %%% Create the 'run' shell file
  fid = fopen(fullfile(resultspath,'run.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,runcommands);
  fclose(fid);
  
  %%% Copy other files across
  resultslist = dir(['./' select_DEFAULTS '/results/']);
  for n=1:1:length(resultslist)
    %%% Ignore hidden files and run script template
    if ((resultslist(n).name(1) == '.') || strcmp(resultslist(n).name,'run_mitgcm'))
      continue;
    end    
    copyfile(fullfile(['./' select_DEFAULTS '/results/'],resultslist(n).name),fullfile(resultspath,resultslist(n).name));
  end    


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILE FOR COMPILING AND RUNNING MITGCM %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% A file to build and run MITgcm, for user convenience if they're
  %%% confident that everything is set up correctly

  %%% Commands to link input and build folders to results folder
  commands = [...
    'cd ./',builddir,'/',lf, ...
    'sh build.sh',lf, ...
    'cd ../',resultsdir,'/ ',lf, ...
    'sh run.sh',lf ];

  %%% Create the 'run' shell file
  fid = fopen(fullfile(exppath,'build_and_run.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,commands);
  fclose(fid);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE SHELL FILES FOR UPLOADING AND DOWNLOADING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Upload command
  upcommand = [...
    'rsync -av --update ',...
    '../',exp_name,' ', ...
    username,'@',clustername,':',clusterdir];
  fid = fopen(fullfile(exppath,'upload_to_cluster.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,upcommand);
  fclose(fid);
  
  %%% Download command
  downcommand = [...
    'rsync -av --update ', ...
    username,'@',clustername,':', ...
    fullfile(clusterdir,exp_name,resultsdir),'/* ', ...
    './results/ \n'];  
  fid = fopen(fullfile(exppath,'download_from_cluster.sh'),'w');
  if (fid == -1)
    error('Could not open run script for writing');
  end
  fprintf(fid,downcommand);  
  fclose(fid);

  %%% Save output data 
    save([exppath '/setParams.mat']);
end


