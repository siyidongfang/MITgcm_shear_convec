%%%
%%% setParams.m
%%%
%%% Sets basic MITgcm parameters plus parameters for included packages, and
%%% writes out the appropriate input files.
%%%

function [nTimeSteps,h,tNorth,sNorth,rho_north,N]...
    = setParams(exp_name,inputpath,codepath,imgpath,listterm,Nx,Ny,Nr,Atide,randtopog_height,randtopog_length,run_type,Shear)

  FigureIsVisible = true;
  addpath ../utils/;
  addpath ../newexp_utils/;
  addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/;
  addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/;
  addpath /Users/ysi/Software/gsw_matlab_v3_06_11/;

  %%%%%%%%%%%%%%%%%%
  %%%%% SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%      
  
  %%% If set true, plots of prescribed variables will be shown
  showplots = true; 
  fontsize = 16;
  fignum = 1;
  
  %%% Data format parameters
  ieee='b';
  prec='real*8';
  realdigits = 8;
  realfmt=['%.',num2str(realdigits),'e'];
  
  %%% Get parameter type definitions
  paramTypes;     

  %%% To store parameter names and values
  parm01 = parmlist;
  parm02 = parmlist;
  parm03 = parmlist;
  parm04 = parmlist;
  parm05 = parmlist;
  PARM={parm01,parm02,parm03,parm04,parm05}; 
  
  %%% Seconds in one hour
  t1min = 60;
  %%% Seconds in one hour
  t1hour = 60*t1min;
  %%% hours in one day
  t1day = 24*t1hour;
  %%% Seconds in 1 year
  t1year = 365*t1day;  
  %%% Metres in one kilometre
  m1km = 1000; 
  %%% Pascals in 1 decibar
  Pa1dbar = 1e4;
      
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% FIXED PARAMETER VALUES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  simTime = 20*t1day;
   % simTime = 1000;
  nIter0 = 0;
  % if(run_type=='init')
  %     simTime = 1*t1day;
  %     nIter0 = 0; %%% Initial iteration 
  %     dt0 = 2;
  % elseif(run_type=='spin')
  %     simTime = 15*t1day;
  %     dt0 = 5;
  %     % nIter0 = 1*t1day/dt0;
  %     nIter0 = 354240;
  % elseif(run_type=='prod')
  %     % simTime = 5*t1day+1*t1hour;
  %     % dt0 = 5;
  %     % nIter0 = 25*t1day/dt0;
  %     simTime = 15*t1day;
  %     nIter0 = 0;
  % end
  
  
  Ly = 3*m1km;
  Lx = 3*m1km; 

  g = 9.81; %%% Gravity
  Omega = 2*pi*366/365/86400;
  Rp = 6400*m1km; %%% Planetary radius

  %%%%%%%%%% The following code calculate the Coriolis parameter:  
  %%%%%%%%%% We know the Coriolis parameter of the observation and the frequency of the M2 tide omega_M2.
  %%%%%%%%%% In the idealized simulations we use an idealized tidal period, which is different from omega_M2
  %%%%%%%%%% We want to keep the slope of internal wave characteristics r_iw (Lamb 2013, Annual Reviews, Eq.3)
  %%%%%%%%%% of the simulation the same as the observation. To keep r_iw the same, so we change the Coriolis parameter for the simulations. 
  omega_M2 = 2*pi/44712;
  omega_idealized = 2*pi/43200;
  latN = 54.2;
  % f_obs = 4*pi/86164*sind(54.2); 
  % load CTD.mat
  % N_const = sqrt(mean(N2_mean15(750:end)));
  % % N_const = 0.62*2.4/1e3; %%% Median buoyancy frequency from Hans van Haren, et al. (preprint)
  % 
  % N1 = N_const;N2 = N_const;
  % f1 = f_obs;
  % omega1 = omega_M2;omega2 = omega_idealized;
  % f2_squared = omega2^2 - (N2^2-omega2^2)*(omega1^2-f1^2)/(N1^2-omega1^2);

  % f0 = sqrt(f2_squared)

  %%%%%%%%%% end of calculating the Coriolis parameter   
  % f0 = 0.53e-4; %%% Coriolis parameter         %-- from Xiaozhou, reference latitude = asind(0.53/10000*86400/4/pi) = 21.37 degree N
  % beta = 1e-11; %%% Beta parameter    
  f0 = 1.18e-4;
  % f0 = 0;
  beta = 0;                                    %-- from Xiaozhou
  rhoConst = 999.8; %%% Reference density       %-- from Xiaozhou, MITgcm default value 999.8

  nonHydrostatic = true; 

  varyingtidalphase = false; % Set true to include zonally (along-slope) varying tidal phase 
  useLAYERS = false;
  % if(run_type=='prod')
  %     useLAYERS = true;
  % end
  useEXF = false;

  % useXRUANcode = false;
  useOBCS = false;
  useRBCS = false;
  %%% OBCS package options
  % if(useXRUANcode)
  %     useOBCS = false;
  %     useRBCS = false;
  % else
  %     useOBCS = true;
  %     useRBCS = false;
  % end
  
  if(useOBCS)
      useOBCStides = false;
      if(Atide~=0)
         useOBCStides = true;
      end

    useobcsNorth = true;
    useobcsSouth = true;
    useOrlanskiNorth = false;
    useOrlanskiSouth = false;
  end

  % Zonal boundary condition
  usePeriodic = true;       %%% Periodic boundary condition for the zonal boundaries
  if(usePeriodic)
      useobcsEast = false;
      useobcsWest = false; 
  end

  %%% RBCS

  %%% Flag for barotropic mode
  isBarotropic = Nr == 1;
  
  %%% PARM01
  %%% momentum scheme
  % parm01.addParm('vectorInvariantMomentum',true,PARM_BOOL); %%% test20231027
  % parm01.addParm('implicSurfPress',0.6,PARM_REAL); %%% test20231027
  % parm01.addParm('implicDiv2DFlow',0.6,PARM_REAL); %%% test20231027

  %------ ysi's configuration
  % viscAh = 1e-4; %%% Horizontal viscosity 
  % viscAr = 2e-4; %%% Vertical viscosity 
  % diffKhT = 1e-5; %%% Horizontal temp diffusion
  % diffKhS = 1e-5; %%% Horizontal salt diffusion
  % diffKrT = 2e-5; %%% Vertical temp diffusion   
  % diffKrS = 2e-5; %%% Vertical salt diffusion 

  %------ xruan's viscosity and diffusivity
  lfac = 2; 
  viscAh = 1e-4*lfac; %%% Horizontal viscosity         %-- from Xiaozhou
  viscAr = 2e-4*lfac; %%% Vertical viscosity           %-- from Xiaozhou
  diffKhT = 1e-4*lfac; %%% Horizontal temp diffusion   %-- from Xiaozhou
  diffKhS = 1e-4*lfac; %%% Horizontal salt diffusion   %-- from Xiaozhou
  diffKrT = 2e-4*lfac; %%% Vertical temp diffusion     %-- from Xiaozhou
  diffKrS = 2e-4*lfac; %%% Vertical salt diffusion     %-- from Xiaozhou
  %------ xruan's viscosity and diffusivity

  viscA4 = 0; %%% Biharmonic viscosity
  viscAhGrid = 0; %%% Grid-dependent viscosity
  viscA4Grid = 0;   
  viscC4smag = 4;  
  diffK4Tgrid = 0; 
  diffK4Sgrid = 0; 
  parm01.addParm('viscA4',viscA4,PARM_REAL);
  parm01.addParm('viscA4Grid',viscA4Grid,PARM_REAL);
  parm01.addParm('viscAhGrid',viscAhGrid,PARM_REAL);
  parm01.addParm('viscA4GridMax',0.5,PARM_REAL);
  parm01.addParm('viscAhGridMax',1,PARM_REAL);
  parm01.addParm('useAreaViscLength',false,PARM_BOOL);
  parm01.addParm('useFullLeith',true,PARM_BOOL);
  parm01.addParm('viscC4smag',viscC4smag,PARM_REAL);  
  % parm01.addParm('useSmag3D',true,PARM_BOOL);     %%% 2023-08-10
  % parm01.addParm('smag3D_coeff',1.0e-2,PARM_REAL);%%% 2023-08-10, MITgcm default value: 1.0e-2
  parm01.addParm('viscC4leith',0,PARM_REAL);
  parm01.addParm('viscC4leithD',0,PARM_REAL);  
  parm01.addParm('viscC2leith',0,PARM_REAL);
  parm01.addParm('viscC2leithD',0,PARM_REAL);  
  %------ ysi's configuration


  parm01.addParm('viscAr',viscAr,PARM_REAL);      
  parm01.addParm('viscAh',viscAh,PARM_REAL);
  parm01.addParm('diffKrT',diffKrT,PARM_REAL);
  parm01.addParm('diffKhT',diffKhT,PARM_REAL);
  parm01.addParm('diffKrS',diffKrS,PARM_REAL);
  parm01.addParm('diffKhS',diffKhS,PARM_REAL);

  %%% diffusivity
  % parm01.addParm('tempAdvScheme',80,PARM_INT);
  % parm01.addParm('saltAdvScheme',80,PARM_INT);
  parm01.addParm('saltAdvScheme',7,PARM_INT);
  parm01.addParm('tempAdvScheme',7,PARM_INT); %--xiaozhou: 33
  parm01.addParm('tempStepping',true,PARM_BOOL); 
  parm01.addParm('saltStepping',true,PARM_BOOL);
  parm01.addParm('staggerTimeStep',true,PARM_BOOL);
  %%% equation of state
  % parm01.addParm('eosType','MDJWF',PARM_STR); 
  parm01.addParm('eosType','LINEAR',PARM_STR); 
  tAlpha = 2e-4;
  sBeta = 0;
  tRef = 0;
  sRef = 35;
  parm01.addParm('tAlpha',tAlpha,PARM_REAL); %-- default
  parm01.addParm('sBeta',sBeta,PARM_REAL);   %-- from Xiaozhou
  parm01.addParm('tRef',tRef,PARM_REAL);     %-- from Xiaozhou
  parm01.addParm('sRef',sRef,PARM_REAL);     %-- from Xiaozhou
  %%% boundary conditions
  parm01.addParm('no_slip_sides',false,PARM_BOOL);
  parm01.addParm('no_slip_bottom',false,PARM_BOOL);
  parm01.addParm('bottomDragLinear',0,PARM_REAL);  %-- from Xiaozhou
  % parm01.addParm('bottomDragQuadratic',2.5e-3,PARM_REAL);  %-- from Xiaozhou
  parm01.addParm('bottomDragQuadratic',0,PARM_REAL); 
  %%% physical parameters
  parm01.addParm('f0',f0,PARM_REAL);
  parm01.addParm('beta',beta,PARM_REAL);
  parm01.addParm('gravity',g,PARM_REAL);
  %%% full Coriolis force parameters
  % parm01.addParm('quasiHydrostatic',false,PARM_BOOL);  %-- from xiaozhou
  % parm01.addParm('fPrime',0,PARM_REAL);  %-- from xiaozhou
  % parm01.addParm('rhoConst',rhoConst,PARM_REAL); %-- from xiaozhou
  %%% implicit diffusion and convective adjustment  
  parm01.addParm('ivdc_kappa',0,PARM_REAL); %%% reference value, 1.0 implicit vertical diffusivity for convection (m2/s)
  parm01.addParm('implicitDiffusion',true,PARM_BOOL);
  parm01.addParm('implicitViscosity',true,PARM_BOOL);
  %%% exact volume conservation
  % parm01.addParm('exactConserv',true,PARM_BOOL);
  parm01.addParm('exactConserv',false,PARM_BOOL); %-- from xiaozhou test20231027
  %%% C-V scheme for Coriolis term
  % parm01.addParm('useCDscheme',false,PARM_BOOL);%-- from xiaozhou
  %%% partial cells for smooth topography
  if (isBarotropic)
    parm01.addParm('hFacMin',0,PARM_REAL);  
  else
    % parm01.addParm('hFacMin',0.1,PARM_REAL);  
    %%% TO DO: set the minimum fraction size of a cell to be 1/10 of the vertical
    %%% grid spacing near the seafloor.
    parm01.addParm('hFacMin',0.3,PARM_REAL);  %-- from Xiaozhou
  end
  %%% file IO stuff
  parm01.addParm('readBinaryPrec',64,PARM_INT);
  parm01.addParm('writeBinaryPrec',64,PARM_INT);
  parm01.addParm('useSingleCpuIO',true,PARM_BOOL);
  parm01.addParm('debugLevel',-1,PARM_INT);
  %%% Wet-point method at boundaries - may improve boundary stability
  parm01.addParm('useJamartWetPoints',true,PARM_BOOL); 
  parm01.addParm('useJamartMomAdv',true,PARM_BOOL); 
  % parm01.addParm('rhoConst',1000,PARM_REAL);
  parm01.addParm('useRealFreshWaterFlux',false,PARM_BOOL); 
  %%% useRealFreshWaterFlux: use true E-P-R freshwater flux (changes free
  %%% surface/sea level) on/off flag, default: false
  parm01.addParm('nonHydrostatic',nonHydrostatic,PARM_BOOL);  
  parm01.addParm('momTidalForcing',false,PARM_BOOL);  




  %%% PARM02
  % parm02.addParm('useSRCGSolver',true,PARM_BOOL);  %-- from Xiaozhou
  % parm02.addParm('cg2dMaxIters',1000,PARM_INT);  
  parm02.addParm('cg2dMaxIters',10000,PARM_INT);          %-- from Xiaozhou
  % parm02.addParm('cg2dTargetResidual',1e-12,PARM_REAL);
  parm02.addParm('cg2dTargetResidual',1e-14,PARM_REAL);   %-- from Xiaozhou
  parm02.addParm('cg3dMaxIters',400,PARM_INT);            %-- from Xiaozhou
  parm02.addParm('cg3dTargetResidual',1e-14,PARM_REAL);   %-- from Xiaozhou

  %%% PARM03
  % parm03.addParm('alph_AB',1/2,PARM_REAL); %-- from Xiaozhou
  % parm03.addParm('beta_AB',5/12,PARM_REAL);%-- from Xiaozhou
  parm03.addParm('forcing_In_AB',false,PARM_BOOL); 
      % This flag makes to model do a  separate (Eulerian?) time step 
      % for the tendencies due to surface forcing. This is sometime 
      % favorable for stability reasons (and some package such as 
      % seaice work only with this).
  parm03.addParm('nIter0',nIter0,PARM_INT);
  parm03.addParm('abEps',0.1,PARM_REAL);
  parm03.addParm('chkptFreq',60*t1hour,PARM_REAL); % rolling 
  parm03.addParm('pChkptFreq',60*t1hour,PARM_REAL); % permanent
  parm03.addParm('taveFreq',0,PARM_REAL); % it only works properly, if taveFreq is a multiple of the time step deltaT (or deltaTclock).
  parm03.addParm('dumpFreq',60*t1hour,PARM_REAL); % interval to write model state/snapshot data (s)
  parm03.addParm('monitorFreq',60*t1hour,PARM_REAL); % interval to write monitor output (s)
  parm03.addParm('dumpInitAndLast',true,PARM_BOOL);
  parm03.addParm('pickupStrictlyMatch',false,PARM_BOOL); 
  parm03.addParm('cAdjFreq',0,PARM_REAL); %%% set to -1, frequency of convective adj. scheme == deltaT

  % %%% Periodic Forcing
  % parm03.addParm('periodicExternalForcing',true,PARM_BOOL); 


  %%% PARM04
  parm04.addParm('usingCartesianGrid',true,PARM_BOOL);
%   parm04.addParm('usingCurvilinearGrid',true,PARM_BOOL);
  parm04.addParm('usingSphericalPolarGrid',false,PARM_BOOL);    
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% GRID SPACING %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%    

  %%% Zonal grid
  dy = Ly/Ny;  
  yy = (1:Ny)*dy;
  
  % %%% Uniform meridional grid   
  % dy = (Ly/Ny)*ones(1,Ny);  

  % dxConst = 100;
  % dxSponge = 100;
  % % Thickness of horizontal sponge layers in gridpoints  
  dxSponge = Lx/Nx;
  spongeThicknessDim = 1*m1km;
  spongeThickness = round(spongeThicknessDim/dxSponge);
  % 
  % dx = [dxSponge*ones(1,spongeThickness) ...
  %     dxConst*ones(1,Nx-2*spongeThickness) ...
  %     dxSponge*ones(1,spongeThickness)]; 

  dx = Lx/Nx*ones(1,Nx);
  xx = cumsum((dx + [0 dx(1:end-1)])/2);

  %%% Plotting mesh
  [Y,X] = meshgrid(yy,xx);
  
  % % dz = [1*ones(1,250) [1:3/99:4] 4*ones(1,100)];
  % dz = [1*ones(1,300) [1:3/99:4] 4*ones(1,100)];
  % dz = flipud(dz')';

  dz_const = 3;
  dz = dz_const*ones(1,Nr);

  % % %%% Varied dz with depth  %  -- from Xiaozhou
  % % % Hsurface = 1002;
  % % % Ntop = 120;
  % % % dz = dz_const.*ones(1,Nr);
  % % % dz(Nr-Ntop + 1:Nr) = dz(Nr - Ntop) * 1.015.^(1:Ntop);
  % % % sum_dz_sponge = sum(dz(Nr-Ntop + 1:Nr));
  % % % dz(Nr-Ntop + 1:Nr) = dz(Nr-Ntop + 1:Nr).*Hsurface/sum_dz_sponge;
  % % % dz = flipud(dz')';
 
  zz = -cumsum((dz+[0 dz(1:end-1)])/2);


  if (length(zz) ~= Nr)
    error('Vertical grid size does not match vertical array dimension!');
  end


  %%%%%% Flat bottom -- start
  % Hmax = 950;
  Hmax = 1500;
  h = -Hmax*ones(Nx,Ny);
  %%%%%% Flat bottom -- end

  % %%%%%% Add walls to the flat bottom -- start
  % h(1) = 0;
  % h(end) = 0;
  % %%%%%% Add walls to the flat bottom -- end


% % %%%%%% Single slope -- start
% %     Hdeep = 2000;
% %     Hshallow = 1400;
% % 
% %     delH = Hdeep-Hshallow; %%% elevation of the slope
% %     theta_slope = 4; %%% 4 degree slope
% %     delL = delH/tand(theta_slope); %%% length of the slope
% % 
% %     Lflat = (Ly-delL)/2;
% % 
% %     YUpslopeStart = Lflat;
% %     YUPslopeEnd = Lflat + delL;
% % 
% %     hh = [Hdeep Hdeep Hshallow Hshallow ];
% %     ll = [yy(1) YUpslopeStart YUPslopeEnd yy(end)];
% %     llf = yy;
% %     h = -interp1(ll,hh,llf,'Linear');
% %     h = smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(h))))))))))))))))))))';
% % %%%%%% Single slope -- end


% %%%%%% Double slopes -- start
%     Hmax = 2000;
%     Hshallow = 1000;
%     Lmiddle = 6*m1km;
% 
%     delH = Hmax-Hshallow; %%% elevation of the slope
%     theta_slope = 4; %%% 4 degree slope
%     delL = delH/tand(theta_slope); %%% length of the slope
% 
%     Lflat = (Lx-delL*2-Lmiddle)/2;
% 
%     YUpslopeStart = Lflat;
%     YUPslopeEnd = Lflat + delL;
%     YDownslopeStart = YUPslopeEnd + Lmiddle;
%     YDownslopeEnd = Lx-Lflat;
% 
%     hh = [Hmax Hmax Hshallow Hshallow Hmax Hmax];
%     ll = [xx(1) YUpslopeStart YUPslopeEnd YDownslopeStart YDownslopeEnd xx(end)];
%     llf = xx;
%     h = -interp1(ll,hh,llf,'Linear');
%     h = smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(smooth(h))))))))))';
% %%%%%% Double slopes -- end


% % %%%%%% Sinusoidal topography -- start
%       ky = 1;
%       Mean_h = -1500*ones(1,Ny);
%       h = zeros(1,Ny);
%       % Ly_tmp = Ly - dy(1);
%       % Ly_flat = 1*m1km;
%       % Ly_tmp = Ly - 2*Ly_flat;
%       % Ly_tmp = 7500; %%% defalt: Ly_tmp = 7500;
%       Hsill = 400; %-- height of the sill
%       Wcanyon = Ly-200;
%       % Wcanyon = 4000; %%% Width of the canyon, default: 7500 m; narrow canyon 4000 m
%       Ly_flat = (Ly - Wcanyon)./2;
%       Ny_flat = round(Ly_flat/dy(1));
%       yidx_flat = [1:Ny_flat-1 Ny-Ny_flat+1:Ny];
%       yidx_sill = [Ny_flat: Ny-Ny_flat];
%       for j = yidx_sill
%           h(1,j) = Mean_h(1,j) +  (Hsill/2 * cos((2*pi/(Wcanyon))*ky .* (yy(j)-Ly_flat)));
%       end
%       for j=yidx_flat
%           h(1,j) = h(1,Ny_flat);
%       end
% % %%%%%% Sinusoidal topography -- end


%%%%%% A simple sloping topography -- start
    % load topog1D.mat   
    % 
    % mean_slope1 = (depth9(3)-depth9(5))/(along_canyon(5)-along_canyon(3))/1000;
    % mean_slope2 = (depth9(5)-depth9(9))/(along_canyon(9)-along_canyon(5))/1000;
    % 
    % Hmid = depth9(5);
    % Hdeep = Hmid + Ly/2*mean_slope1;
    % Hshallow = Hmid - Ly/2*mean_slope2;
    % 
    % hh = [Hdeep Hmid Hshallow];
    % ll = [yy(1) Ly/2 yy(end)];
    % llf = yy;
    % h = -interp1(ll,hh,llf,'Linear');
%%%%%% A simple sloping topography -- end

% %%% Tanh shape  -- start
%     Hs = 1420;
%     Ys = Ly/2;
%     Zs = 1600+200;
%     Ws = Ly/4;
%     h = -(Zs+(Hs/2)*tanh(-(Y-Ys)/Ws));  
% %%% Tanh shape  -- end
% 
% %%%%%% Real topography of the center of the Rockall canyon -- start
%     load topog1D_new.mat   
%     hh = depthn;
%     ll = along_canyonn*1000 - (along_canyonn(end)*1000-Ly);
%     llf = yy;
%     h = -interp1(ll,hh,llf,'Linear');
% 
%     [mm yidx_slope] = min(abs(h+Hs));
%     delH = h2(yidx_slope) - h(yidx_slope);
%     h(yidx_slope:end) = h2(yidx_slope:end)-delH;
%     h(yidx_slope-5:yidx_slope+5) = smooth(smooth(smooth(smooth(smooth(smooth(h(yidx_slope-5:yidx_slope+5)))))))';
% %%%%%% Real topography  -- end


  if(randtopog_length~=0 && randtopog_height~=0)
      %---- Add random bumps to the topography
      h_rand = genRandField_y(randtopog_length,[],randtopog_height,Ny,Ly);
      h_rand = h_rand - min(min(h_rand));
      % h_rand = h_rand .* (-h(2:end)/H); %%% Scale by topographic height so that bumps are most pronounced in deep areas
      % h(2:end) = h(2:end) + h_rand;
      H = sum(dz);
      h_rand = h_rand .* (-h/H); %%% Scale by topographic height so that bumps are most pronounced in deep areas
      % h_rand(1) = 0;
      % h_rand(end) = 0;
      aaa = (h_rand(end)+h_rand(1))/2;
      h_rand(1) = aaa;
      h_rand(end) = aaa;
      h = h + h_rand;
      %----------------
  end

  % %%% Adjust zz according to the deepest topography
  % if (abs(zz(end))< abs(min(h)))   % if (abs(zz(end))< abs(min(h)-dz_const/2))
  %    frac_h2zz = (min(h)-dz_const)./zz(end)
  %    dz = frac_h2zz.*dz;
  %    zz = -cumsum((dz+[0 dz(1:end-1)])/2);
  % end



  %%% Store grid spacings
  parm04.addParm('delX',dx,PARM_REALS);
  parm04.addParm('delY',dy*ones(1,Ny),PARM_REALS);
  parm04.addParm('delR',dz,PARM_REALS);      
  
  % %%% Don't allow partial cell height to fall below min grid spacing
  % if (~isBarotropic)
  %   parm01.addParm('hFacMinDr',min(dz),PARM_REAL);  %-- from Xiaozhou
  % end


  %%% Plot bathymetry
  h_figure=figure(fignum);
  fignum = fignum + 1;
  set(h_figure,'Visible', FigureIsVisible);clf;
  fontsize = 15;
  plot(xx/m1km,h,'LineWidth',2)
  xlabel('Latitude, y (km)')
  ylabel('z (m)')
  title('Bathymetry')
  set(gca,'fontsize',fontsize);grid on;grid minor;
  PLOT = gcf;
  PLOT.Position = [263 149 567 336];
  ylim([-2030 0])
  %%% Save the figure
  savefig([imgpath '/bathymetry.fig']);
  saveas(gcf,[imgpath '/bathymetry.png']);
  
  %%% Save as a parameter
  writeDataset(h,fullfile(inputpath,'bathyFile.bin'),ieee,prec);
  parm05.addParm('bathyFile','bathyFile.bin',PARM_STR); 
  

  
   
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % %%%%% VERTICAL DIFFUSIVITY %%%%%
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % 
  % kappa_max = 5e-3;
  % %%% set the diffusivity here
  % diffKr = diffKrT*ones(Nx,Ny,Nr);  
  % kap_deep = kappa_max; %%% Max diffisivity in deep ocean  
  % H_kap = 300; %%% e-folding scale for mixing decrease with h.a.b.
  % for i=1:Nx
  %   for j=1:Ny      
  %     kap_profile = diffKrT + kap_deep*min(exp(-(zz-h(i,j))/H_kap),1);      
  %     diffKr(i,j,:) = reshape(kap_profile,[1 1 Nr]); 
  %   end
  % end
  % 
  % %%% Plot the surface heat flux
  % if (showplots)
  %   h_figure=figure(fignum);
  %   fignum = fignum + 1;
  %   set(h_figure,'Visible', FigureIsVisible);clf;
  %   set(gcf,'Color','w')
  %   % semilogx(squeeze(diffKr(1,round(Ny/2),:)),-zz); axis ij;grid on;grid minor
  %   plot(squeeze(diffKr(1,round(Ny/2),:)),-zz); axis ij;grid on;grid minor
  %   % title('log(\kappa)');
  %   set(gca,'Fontsize',fontsize)
  %   title('\kappa','Fontsize',fontsize+3)
  %   ylabel('HAB (m)')
  %   xlabel('(m^2/s)')
  % end  
  % 
  % %%% Save as a parameter
  % writeDataset(diffKr,fullfile(inputpath,'diffKrFile.bin'),ieee,prec);
  % parm05.addParm('diffKrFile','diffKrFile.bin',PARM_STR);
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% NORTHERN TEMPERATURE/SALINITY PROFILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

  % Hz = sum (dz);
  % tNorth = (0.62*2.4)^2*(1e-3)^2 *(zz+Hz) /9.81/2e-4;
  % N2const = (1e-3)^2*1e-4;
  % N2const = 1e-10;
  % tNorth = N2const *(zz+Hz) /9.81/2e-4;
  tNorth = 0.*ones(1,Nr);
  sNorth = 0.*ones(1,Nr);

  % tNorth = interp1(pp15',pt_mean15',-zz,'Linear');
  % sNorth = interp1(pp15',psal_mean15',-zz,'Linear');

  %%% Calculate the relaxation density at Northern boundary using GSW toolbox
  ref_pres_surf = 0; 
  lonN = -11.9;
  SA_north = gsw_SA_from_SP(sNorth,ref_pres_surf,lonN,latN);  
  CT_north = gsw_CT_from_pt(SA_north,tNorth); 

  tSouth = tNorth; %%% Southern boundary conditions are the same as the northern boundary
  sSouth = sNorth;

  if(sNorth==0)
      parm05.addParm('checkIniSalt',false,PARM_BOOL);
  end
  if(tNorth==0)
      parm05.addParm('checkIniTemp',false,PARM_BOOL);
  end

  % tNorth = -1.5.*ones(1,Nr);
  % sNorth = 33.7:0.1/(Nr-1):33.8;  %%% weak stratification
  % sNorth = 34.*ones(1,Nr);   %%% zero stratification
  % sNorth = 33:1.5/(Nr-1):34.5; %%% strong stratification

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Calculate density and make plots %%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%% Linear EOS
  rho_north = rhoConst.*(1-(tNorth-tRef)*tAlpha);
    % rho_north = gsw_rho(SA_north,CT_north,ref_pres_surf); 
    % rho_north = densmdjwf(sNorth,tNorth,ref_pres_surf)';
    %%% Plot the relaxation density
  if (showplots)
    h_figure=figure(fignum);
    fignum = fignum + 1;
    set(h_figure,'Visible', FigureIsVisible);clf;
    plot(rho_north,-zz,'LineWidth',1.5); axis ij;
    xlabel('\rho_r_e_f (\circC)');
    ylabel('Depth (m)');
    title('Relaxation density (P_{ref} = 0)');
    set(gca,'fontsize',fontsize);
    PLOT = gcf;
    PLOT.Position = [644 148 380 562];  
    grid on;grid minor;
    %%% Save the figure
    savefig([imgpath '/RelaxationDensity_surf.fig']);
    saveas(gcf,[imgpath '/RelaxationDensity_surf.png']);
  end
   
  
  %%% Plot the relaxation temperature
  if (showplots)
    h_figure=figure(fignum);
    fignum = fignum + 1;
    set(h_figure,'Visible', FigureIsVisible);clf;
    plot(tNorth,-zz,'LineWidth',1.5); axis ij;
    xlabel('\theta_r_e_f (\circC)');
    ylabel('Depth (m)');
    title('Relaxation temperature');
    set(gca,'fontsize',fontsize);
    PLOT = gcf;
    PLOT.Position = [644 148 380 562];  
    grid on;grid minor;
    %%% Save the figure
    savefig([imgpath '/RelaxationT.fig']);
    saveas(gcf,[imgpath '/RelaxationT.png']);
  end
    
    
  %%% Plot the relaxation salinity
  if (showplots)
    h_figure=figure(fignum);
    fignum = fignum + 1;
    set(h_figure,'Visible', FigureIsVisible);clf;
    plot(sNorth,-zz,'LineWidth',1.5);axis ij;
    xlabel('S_r_e_f (psu)');
    ylabel('Depth (m)');
%     ylabel('z','Rotation',0);
    title('Relaxation salinity');
    set(gca,'fontsize',fontsize);
    PLOT = gcf;
    PLOT.Position = [644 148 380 562]; 
    grid on;grid minor;
    %%% Save the figure
    savefig([imgpath '/RelaxationS.fig']);
    saveas(gcf,[imgpath '/RelaxationS.png']);
  end
    
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFORMATION RADIUS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Check Brunt-Vaisala frequency using linear EOS
    % [N2_north, pp_mid_north] = gsw_Nsquared(SA_north,CT_north,-zz,latN);
    N2_north = -g/rhoConst.*(rho_north(1:end-1)-rho_north(2:end))./(zz(1:end-1)-zz(2:end));
    pp_mid_north = 0.5*(-zz(1:end-1)+(-zz(2:end))); %%% Mid-depth where the buoyancy frequency is defined
    dzData = zz(1:end-1)-zz(2:end);
    
    % find(N2_north<0) 

    %%% Calculate internal wave speed and first Rossby radius of deformation
    N = sqrt(N2_north);
    Cig = zeros(size(xx));
    for j=1:Nx    
      for k=1:length(dzData)
        if (zz(k) > h(j,:)) 
            % k
            % min(dzData(k),zz(k)-h(1,j))
          Cig(j) = Cig(j) + N(k)*min(dzData(k),zz(k)-h(j,:));
           % Cig(j) = Cig(j) + N(k)*(zz(k)-h(1,j));
        end
      end
    end
    Rd = Cig./(pi*abs(f0+beta*Y(1,:)));

    if (showplots)
      h_figure=figure(fignum);
      fignum = fignum + 1;
      set(h_figure,'Visible', FigureIsVisible);clf;
      semilogx(N2_north,pp_mid_north,'LineWidth',1.5);axis ij;
      xlabel('N^2 (s^-^2)');
      ylabel('Depth (m)');
%       ylabel('z (km)','Rotation',0);
      % xlim([0.999 1.001]*1e-6)
      % xlim([1 2]*1e-3)
      title('N^2');
      set(gca,'fontsize',fontsize);
      PLOT = gcf;
      PLOT.Position = [644 148 380 562];  
      grid on;grid minor;
      %%% Save the figure
      savefig([imgpath '/BuoyancyFrequency.fig']);
      saveas(gcf,[imgpath '/BuoyancyFrequency.png']);
    end
    

    if (showplots)
      h_figure=figure(fignum);
      fignum = fignum + 1;
      set(h_figure,'Visible', FigureIsVisible);clf;
      plot(xx/1000,Rd/1000,'LineWidth',1.5);
      xlabel('Offshore distance (km)');
      ylabel('R_d (km)');
      title('First baroclinic Rossby deformation radius');
      set(gca,'fontsize',fontsize-1);
      PLOT = gcf;
      PLOT.Position = [263 149 567 336];
      grid on;grid minor;
      %%% Save the figure
      savefig([imgpath '/R_d.fig']);
      saveas(gcf,[imgpath '/R_d.png']);
    end
    
  

  %%% Calculate the slope of internal wave characteristics r_iw (Lamb 2013, Annual Reviews, Eq. 3)
  %%% Calculate the topographic slope s_topog
  %%% Use the mean buoyancy frequency N 100m above the topography to calculate r_iw
  %%% Compare r_iw with s_topog


  % % omega_tides = omega_idealized;
  % % % N_mean = mean(N);
  % % % r_iw = sqrt((omega_tides^2-f0^2)/(N_mean^2-omega_tides^2))
  % % r_iw = NaN*zeros(1,Ny-1);
  % % h_mid = 0.5*(h(1:end-1)+h(2:end));
  % % for i=1:Ny-1
  % %     clear zidx_100 zidx1
  % %     [a,zidx1] = min(abs(-h_mid(i)-p_mid15));
  % %     zidx_100m = zidx1-50:zidx1;
  % %     N_mean(i) = sqrt(mean(N2_mean15(zidx_100m)));
  % %     r_iw(i) = sqrt((omega_tides^2-f0^2)/(N_mean(i)^2-omega_tides^2));
  % % end
  % % 
  % % s_topog = diff(h)./dy(1);
  % % yy_mid = 0.5*(yy(2:end)+yy(1:end-1));
  % % 
  % % h_subcritical = NaN.*zeros(1,Ny);
  % % h_supercritical = NaN.*zeros(1,Ny);
  % % % h_subcritical(find(abs(s_topog)<r_iw))=h(find(abs(s_topog)<r_iw));
  % % % h_supercritical(find(abs(s_topog)>=r_iw))=h(find(abs(s_topog)>=r_iw));
  % % 
  % % for i=1:Ny-1
  % %     if(abs(s_topog(i))<r_iw(i))
  % %         h_subcritical(i) = h(i);
  % %     else
  % %         h_supercritical(i) = h(i);
  % %     end
  % % end
  % % 
  % % 
  % %   if (showplots)
  % %     h_figure=figure(fignum);
  % %     fignum = fignum + 1;
  % %     set(h_figure,'Visible', FigureIsVisible);clf;
  % %     plot(yy_mid/1000,s_topog,'LineWidth',1.5);
  % %     hold on;
  % %     % plot(yy_mid/1000,r_iw.*ones(1,length(yy_mid)),'r--','LineWidth',1.5);
  % %     plot(yy_mid/1000,r_iw,'r--','LineWidth',1.5);
  % %     hold off;
  % %     xlabel('Offshore distance (km)');
  % %     ylabel('S_{topo}');
  % %     title('Topographic Slope');
  % %     set(gca,'fontsize',fontsize-1);
  % %     PLOT = gcf;
  % %     PLOT.Position = [263 149 567 336];
  % %     grid on;grid minor;
  % %     %%% Save the figure
  % %     savefig([imgpath '/s_topog.fig']);
  % %     saveas(gcf,[imgpath '/s_topog.png']);
  % % 
  % % 
  % %     %%% Plot bathymetry
  % %     h_figure=figure(fignum);
  % %     fignum = fignum + 1;
  % %     set(h_figure,'Visible', FigureIsVisible);clf;
  % %     plot(yy/m1km,h_subcritical,'LineWidth',2)
  % %     hold on;
  % %     plot(yy/m1km,h_supercritical,'LineWidth',2)
  % %     hold off;
  % %     legend('Subcritical','Supercritical','Position',[0.4083 0.7515 0.2055 0.1086])
  % %     xlabel('Latitude, y (km)')
  % %     ylabel('z (m)')
  % %     title('Bathymetry')
  % %     set(gca,'fontsize',fontsize);grid on;grid minor;
  % %     PLOT = gcf;
  % %     PLOT.Position = [263 149 567 336];
  % %     xlim([0 Ly/1000])
  % %     % ylim([-1800 -1200])
  % %     %%% Save the figure
  % %     savefig([imgpath '/bathymetry.fig']);
  % %     saveas(gcf,[imgpath '/bathymetry.png']);
  % % 
  % %   end
  

  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CALCULATE TIME STEP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  
  %%% These estimates are in no way complete, but they give at least some
  %%% idea of the time step needed to keep things stable. In complicated 
  %%% simulations, preliminary tests may be required to estimate the
  %%% parameters used to calculate these time steps.        
  
  %%% Gravity wave CFL

  %%% Upper bound for absolute horizontal fluid velocity (m/s)
  %%% At the moment this is just an estimate
  Umax = 1
  %%% Max gravity wave speed
  cmax = max(Cig)
  %%% Max gravity wave speed using total ocean depth
  cgmax = Umax + cmax;
  %%% Advective CFL
  deltaT_adv = min([0.5*dx/Umax,0.5*dy/Umax]);
  %%% Gravity wave CFL
  deltaT_gw = min([0.5*dx/cmax,0.5*dy/cmax]);
  %%% CFL time step based on full gravity wave speed
  deltaT_fgw = min([0.5*dx/cgmax,0.5*dy/cgmax]);
    
  %%% Other stability conditions
  
  %%% Inertial CFL time step (Sf0<=0.5)
  deltaT_itl = 0.5/abs(f0);
  %%% Time step constraint based on horizontal diffusion 
  deltaT_Ah = 0.5*min([dx dy])^2/(4*viscAh);    
  %%% Time step constraint based on vertical diffusion
  deltaT_Ar = 0.5*min(dz)^2 / (4*viscAr);  
  %%% Time step constraint based on biharmonic viscosity 
  deltaT_A4 = 0.5*min([dx dy])^4/(32*viscA4);
  %%% Time step constraint based on horizontal diffusion of temp 
  deltaT_KhT = 0.4*min([dx dy])^2/(4*diffKhT);    
  %%% Time step constraint based on vertical diffusion of temp 
  deltaT_KrT = 0.4*min(dz)^2 / (4*diffKrT);
 
  %%% Time step size  
  deltaT = min([deltaT_fgw deltaT_gw deltaT_adv deltaT_itl deltaT_Ah deltaT_Ar deltaT_KhT deltaT_KrT deltaT_A4]);
  deltaT = floor(deltaT) 
  % deltaT = floor(deltaT*2/3) 


  % deltaT = 1
  % if(deltaT<5)
  %     deltaT = 5
  % end
  
  % deltaT = 18



  nTimeSteps = ceil(simTime/deltaT);
  simTimeAct = nTimeSteps*deltaT
  
  %%% Write end time and time step size  
  parm03.addParm('endTime',nIter0*deltaT+simTimeAct,PARM_INT);
  parm03.addParm('deltaT',deltaT,PARM_REAL); 


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% TRACER DIFFUSION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set biharmonic diffusivities as a fraction of the maximum stable
  %%% grid-scale hyperdiffusion
  diffK4T = diffK4Tgrid * max([dx dy])^4 / (32*deltaT);
  diffK4S = diffK4Sgrid * max([dx dy])^4 / (32*deltaT);
  % parm01.addParm('diffK4T',diffK4T,PARM_REAL); %-- from Xiaozhou
  % parm01.addParm('diffK4S',diffK4S,PARM_REAL); %-- from Xiaozhou



  
  %%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL DATA %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%
    
  %%% Random noise amplitude
  tNoise = 1e-20;  
  % tNoise = 0;
  sNoise = 0;

  %---- Add an infinitesimal linear stratification 
  Flinear = zeros(1,Nr);
  % Flinear = tNoise*(flip(-zz)); 

  %---- Add random noise with a certain wavelength to the initial temperature field
  noise_length = 400;
  noise_amp = 1;
  Nx_noise = Lx;
  Nr_noise = Hmax;

  Fnoise = genRandField_xz(noise_length,[],noise_amp,Nx_noise,Nr_noise,Lx,Hmax);
  Fnoise = tNoise*Fnoise/max(max(abs(Fnoise)));
  [zzz1,xxx1] = meshgrid(-1*(1:Hmax),1:Lx);
  [zzz2,xxx2] = meshgrid(zz,xx);

  Fnoise = interp2(zzz1,xxx1,Fnoise,zzz2,xxx2);
  %----------------

  %%% Align initial temp with background
  hydroTh = ones(Nx,Ny,Nr);
  hydroSa = ones(Nx,Ny,Nr);

  for k=1:Nr
    hydroTh(:,:,k) = squeeze(hydroTh(:,:,k))*tNorth(k);
    hydroSa(:,:,k) = squeeze(hydroSa(:,:,k))*sNorth(k);
  end

  % %%%% Initial condition: random noise in buoyancy field
  % hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);
  % hydroSa = hydroSa + sNoise*(2*rand(Nx,Ny,Nr)-1);

  for i=1:Nx
      for j=1:Ny
          for k=1:Nr
              % hydroTh(i,j,k)=hydroTh(i,j,k)+Flinear(k);
              % hydroTh(i,j,k)=hydroTh(i,j,k)+Fnoise(i,k)+Flinear(k);
              hydroTh(i,j,k)=hydroTh(i,j,k)+Fnoise(i,k);
          end
      end
  end


  % %%% Titled isotherms
  % theta_slope = 4; %%% 4 degree slope
  % 
  % % N2const = (1e-3)^2;
  % % tNorth = N2const *(zz+Hz) /9.81/2e-4;
  % 
  % for i=1:Nx
  %     for j=1:Ny
  %         for k=1:Nr
  %             hydroTh(i,j,k) = N2const/9.81/2e-4 *...
  %                 ( cosd(theta_slope)*(zz(k)+Hz) + sind(theta_slope)*xx(i));
  %         end
  %     end
  % end

  h_figure=figure(fignum);
  fignum = fignum + 1;
  set(h_figure,'Visible', FigureIsVisible);clf;set(gcf,'Color','w')
  pcolor(xx/1000,-zz,squeeze(hydroTh)');axis ij;
  % hold on;
  % contour(xx/1000,-zz,squeeze(hydroTh)',[0:0.03:2.5]);axis ij;
  hold on;
  plot(xx/1000,-h);
  shading flat;colormap(redblue);
  if(tNoise ~=0)
      if (max(Flinear)~=0)
        clim([-tNoise tNoise]*10000);
      else
        clim([-tNoise tNoise]);
      end
  end
  % clim([0 1]);
  colorbar;
  xlabel('Distance, y (km)');
  ylabel('Depth (m)');
  title('Initial temperature (^oC)');
  set(gca,'fontsize',fontsize);
  PLOT = gcf;
  PLOT.Position = [80 220 463 342];
  grid on;grid minor;
  %%% Save the figure
  savefig([imgpath '/hydroTh.fig']);
  saveas(gcf,[imgpath '/hydroTh.png']);
  
  %%% Add some random noise
  if (~isBarotropic)
    hydroTh = hydroTh + tNoise*(2*rand(Nx,Ny,Nr)-1);
    hydroSa = hydroSa + sNoise*(2*rand(Nx,Ny,Nr)-1);
  end
  
  %%% Write to data files
  writeDataset(hydroTh,fullfile(inputpath,'hydrogThetaFile.bin'),ieee,prec); 
  parm05.addParm('hydrogThetaFile','hydrogThetaFile.bin',PARM_STR);
  writeDataset(hydroSa,fullfile(inputpath,'hydrogSaltFile.bin'),ieee,prec); 
  parm05.addParm('hydrogSaltFile','hydrogSaltFile.bin',PARM_STR); 


  %%% Restore temperature and velocity shear at the horizontal boundaries
  HoriSpongeIdx = [1:spongeThickness Ny-spongeThickness+1:Ny];
  vrelax = zeros(1,Nr);
  % Shear and Hshear must be changed together with external_forcing.F
  h_shear = 250;
  % Hshear = Hmax-h_shear; 
  % [a Nshear] = min(abs(abs(zz)-Hshear));
  % for k = Nshear:Nr
  %     vrelax(k) = (zz(k)-zz(end)+dz_const/2)*Shear;
  % end
  % for k = 1:Nshear-1
  %     vrelax(k) = vrelax(Nshear);
  % end

  % Nshear_smooth_half = round(15*3/dz_const);
  % Nshear_smooth_half = 50;
  Nshear_smooth_half = 40;
  % Nshear_smooth_half = 0;
  % Nsmooth_span = Nshear_smooth_half*2+1;
  % vrelax = smooth(vrelax,Nsmooth_span);


  %%
  % shearProfile = zeros(1,Nr); 
  for i=1:Nr
       if((zz(i)-zz(Nr))<h_shear) 
           shearProfile(i)=(zz(i)+Hmax)/h_shear;
       else
           shearProfile(i)=1.;
       end 
   end 

  vrelax2 = Shear*h_shear*shearProfile;
  

  %--- smooth the velocity shear
  for kLev = 1:Nr
      if(kLev>Nshear_smooth_half) 
          if((Nr-kLev)>=Nshear_smooth_half) 
              NsmoothStart = kLev-Nshear_smooth_half;
              NsmoothEnd = kLev+Nshear_smooth_half;
          end 
      end 

      if((Nr-kLev)<Nshear_smooth_half) 
          NsmoothStart = kLev-(Nr-kLev);
          NsmoothEnd = Nr;
      end 

      if(kLev<=Nshear_smooth_half) 
          NsmoothStart = 1;
          NsmoothEnd = kLev+(kLev-1);
      end 

      shearRatio(kLev) = 0.;
      Ndivide(kLev) = 0.;

      if(NsmoothEnd>NsmoothStart) 
          for i= NsmoothStart:NsmoothEnd
               shearRatio(kLev) = shearRatio(kLev) + shearProfile(i);
               Ndivide(kLev) = Ndivide(kLev) + 1;
          end 
          shearRatio(kLev) = shearRatio(kLev)/Ndivide(kLev);
      else
          shearRatio(kLev) = shearProfile(kLev);
      end 
  end

  vrelax = Shear*h_shear*shearRatio;

  % vrelax = vrelax2;

  %%% Plot velocity shear
  h_figure=figure(fignum);
  fignum = fignum + 1;
  set(h_figure,'Visible', FigureIsVisible);clf;
  fontsize = 15;
  plot(vrelax2,-zz,'LineWidth',2);hold on;
  plot(vrelax,-zz,'LineWidth',2);axis ij;
  xlabel('v (m/s)')
  ylabel('Depth (m)')
  title('Velocity')
  set(gca,'fontsize',fontsize);grid on;grid minor;
  PLOT = gcf;
  PLOT.Position = [263 149 380 552];
  %%% Save the figure
  savefig([imgpath '/vrelax.fig']);
  saveas(gcf,[imgpath '/vrelax.png']);



  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% RBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  
  if(useRBCS)   
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Add a mask at ocean surface to absorb the upward radiating internal
  %%% waves

  %%% To store parameter names and values
  rbcs_parm01 = parmlist;
  rbcs_parm02 = parmlist;
  RBCS_PARM = {rbcs_parm01,rbcs_parm02};
  
  useRBCuVel = true;
  % useRBCvVel = true;
  tauRelaxU = 1*t1min;
  % tauRelaxV = 10*t1min;
  rbcs_parm01.addParm('useRBCuVel',useRBCuVel,PARM_BOOL);
  % rbcs_parm01.addParm('useRBCvVel',useRBCvVel,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxU',tauRelaxU,PARM_REAL);
  % rbcs_parm01.addParm('tauRelaxV',tauRelaxV,PARM_REAL);

  useRBCtemp = true;
  % useRBCsalt = false;
  tauRelaxT = 1*t1min;
  % tauRelaxS = 1*t1hour;
  rbcs_parm01.addParm('useRBCtemp',useRBCtemp,PARM_BOOL);
  % rbcs_parm01.addParm('useRBCsalt',useRBCsalt,PARM_BOOL);
  rbcs_parm01.addParm('tauRelaxT',tauRelaxT,PARM_REAL);
  % rbcs_parm01.addParm('tauRelaxS',tauRelaxS,PARM_REAL);
  

  % temp_relax = zeros(Nx,Ny,Nr);
  % salt_relax = zeros(Nx,Ny,Nr);
  uvel_relax = zeros(Nx,Ny,Nr);
  % vvel_relax = zeros(Nx,Ny,Nr);

  % % %%% Restore surface values 
  % % for i=1:Nx
  % %     for j=1:Ny
  % %         temp_relax(i,j,1:Nsponge_top) = tNorth(1:Nsponge_top); 
  % %         % salt_relax(i,j,1:Nsponge_top) = sNorth(1:Nsponge_top); 
  % %     end
  % % end

  % uvel_relax = 1e-27*random('Normal',0,1).*ones(Nx,Ny,Nr);
  % vvel_relax = 1e-27*random('Normal',0,1).*ones(Nx,Ny,Nr);


  % for i=1:Nx
  %     for j=HoriSpongeIdx
  %         vvel_relax(i,j,:) = vrelax; 
  %     end
  % end


  msk=zeros(Nx,Ny,Nr);
  %%% Horizontal Mask
  for i = 1:spongeThickness
      msk(i,:,:) =(spongeThickness-i+1)./spongeThickness; 
  end
  for i = Nx-spongeThickness+1:Nx
      msk(i,:,:) =(i-(Ny-spongeThickness))./spongeThickness; 
  end

  % % %%% Mask is 1 at ocean surface, and decreases to zero at the bottom of
  % % %%% the sponge layer.
  % % %%% Mask is zero below the sponge layer, i.e. no relaxation
  % % for k=1:Nsponge_top
  % %     msk(:,:,k) =(Nsponge_top-k+1)./Nsponge_top; 
  % % end

  temp_mask = msk; 
  % salt_mask = msk; 
  uvel_mask = msk; 
  % vvel_mask = msk; 

  temp_relax = hydroTh;

  %%% Save as parameters
  writeDataset(temp_relax,fullfile(inputpath,'sponge_temp.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxTFile','sponge_temp.bin',PARM_STR);
  writeDataset(temp_mask,fullfile(inputpath,'rbcs_temp_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskFile(1)','rbcs_temp_mask.bin',PARM_STR);

  % writeDataset(salt_relax,fullfile(inputpath,'sponge_salt.bin'),ieee,prec); 
  % rbcs_parm01.addParm('relaxSFile','sponge_salt.bin',PARM_STR);
  % writeDataset(salt_mask,fullfile(inputpath,'rbcs_salt_mask.bin'),ieee,prec); 
  % rbcs_parm01.addParm('relaxMaskFile(2)','rbcs_salt_mask.bin',PARM_STR);
  % 
  writeDataset(uvel_relax,fullfile(inputpath,'sponge_uvel.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxUFile','sponge_uvel.bin',PARM_STR);
  writeDataset(uvel_mask,fullfile(inputpath,'rbcs_uvel_mask.bin'),ieee,prec); 
  rbcs_parm01.addParm('relaxMaskUFile','rbcs_uvel_mask.bin',PARM_STR);

  % writeDataset(vvel_relax,fullfile(inputpath,'sponge_vvel.bin'),ieee,prec); 
  % rbcs_parm01.addParm('relaxVFile','sponge_vvel.bin',PARM_STR);
  % writeDataset(vvel_mask,fullfile(inputpath,'rbcs_vvel_mask.bin'),ieee,prec); 
  % rbcs_parm01.addParm('relaxMaskVFile','rbcs_vvel_mask.bin',PARM_STR);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.rbcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %%% Creates the 'data.rbcs' file
  write_data_rbcs(inputpath,RBCS_PARM,listterm,realfmt);
  end





  %%%%%%%%%%% Initial data -- part 2

  % if(useXRUANcode)
  %   omega_tides = 1.454441043328608e-4;
  %   vVelInit = 0.*ones(Nx,Ny,Nr);
  %   uVelInit = 0.1 * 0.53e-4/omega_tides * ones(Nx,Ny,Nr);
  %   Rvalue = tNorth';
  %   writeDataset(Rvalue,fullfile(inputpath,'R.init'),ieee,prec); 
  % end

  %%% Initial velocity
  % vVelInit = 0.*ones(Nx,Ny,Nr);
  % uVelInit = 0.*ones(Nx,Ny,Nr);

  % omega = 1.454441043328608e-4;
  % vVelInit = 0.025*ones(Nx,Ny,Nr);


  for i=1:Nx
      for j=1:Ny
          uVelInit(i,j,:) = vrelax*cos(0); 
      end
  end
    
  for i=1:Nx
      for j=1:Ny
          vVelInit(i,j,:) = vrelax*sin(0); 
      end
  end

  h_figure=figure(fignum);
  fignum = fignum + 1;
  set(h_figure,'Visible', FigureIsVisible);clf;set(gcf,'Color','w')
  pcolor(xx/1000,-zz,squeeze(uVelInit)');axis ij;
  hold on;
  plot(xx/1000,-h);
  shading flat;colormap(redblue);clim([-2 2]*(max(vrelax)+0.0001));colorbar;
  xlabel('Distance, y (km)');
  ylabel('Depth (m)');
  title('Initial velocity u (m/s)');
  set(gca,'fontsize',fontsize);
  PLOT = gcf;
  PLOT.Position = [80 220 463 342];
  grid on;grid minor;
  %%% Save the figure
  savefig([imgpath '/uVelInit.fig']);
  saveas(gcf,[imgpath '/uVelInit.png']);

  writeDataset(uVelInit,fullfile(inputpath,'uVelInitFile.bin'),ieee,prec); 
  parm05.addParm('uVelInitFile','uVelInitFile.bin',PARM_STR);
  writeDataset(vVelInit,fullfile(inputpath,'vVelInitFile.bin'),ieee,prec); 
  parm05.addParm('vVelInitFile','vVelInitFile.bin',PARM_STR);
  
  %%% High-resolution runs must be restarted from the end of a
  %%% low-resolution run via doubleRes, which creates the initialization 
  %%% files indicated here
  % if (run_type=='prod')
  %   parm05.addParm('uVelInitFile','uVelInitFile.bin',PARM_STR);  
  %   parm05.addParm('vVelInitFile','vVelInitFile.bin',PARM_STR);  
  %   parm05.addParm('pSurfInitFile','pSurfInitFile.bin',PARM_STR);  %initial free surface position
  % end
   

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Creates the 'data' file
  write_data(inputpath,PARM,listterm,realfmt);
 

  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
 if (useLAYERS)

  %%% To store parameter names and values
  layers_parm01 = parmlist;
  LAYERS_PARM = {layers_parm01};
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% LAYERS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Define parameters for layers package %%%
  

  %%% Number of fields for which to calculate layer fluxes
  % The number of possible layers coordinates (max number of tracer fields used for layer averaging)
  layers_maxNum = 2;
  % layers_maxNum = 1;

%   %%% Specify potential density
  layers_name = char('RHO','TH'); 
  % layers_name = char('TH'); 
%     layers_name = char('RHO'); 

  % layers_bounds = zeros(Nr+2,2);
  layers_bounds = zeros(Nr*2+2,2);

   rhoShift = rhoConst - 1000;
   max_rho = max(rho_north-rhoConst);
   rho_bounds = rhoShift + [min(rho_north-rhoConst)-2 rho_north-rhoConst max(rho_north-rhoConst)+2]; 
   pt_bounds = [min(tNorth)-2 flip(tNorth) max(tNorth)+2];

   [a1 b1] = min(abs(zz+1000));
   c1 = 4;
   d1 = length(1:c1:b1);
   c3 = length(layers_bounds)-d1-1;

   e11 = rho_bounds(b1+1);
   e21 = rho_bounds(end-1)+0.02;

   e12 = pt_bounds(b1+1);
   e22 = pt_bounds(end-1)+0.02;

   layers_bounds(1:d1,1) = rho_bounds(1:c1:b1); %%% use less layers for the top 1100 m
   layers_bounds(1:d1,2) = pt_bounds(1:c1:b1); 
   layers_bounds(d1+1:end-1,1) = e11:(e21-e11)/(c3-1):e21; 
   layers_bounds(d1+1:end-1,2) = e12:(e22-e12)/(c3-1):e22; 
   layers_bounds(end,1) = rho_bounds(end);
   layers_bounds(end,2) = pt_bounds(end);

  %%% Reference level for calculation of potential density
  refDepth = 0*m1km;
  [dz_refDepth idx_refDepth] = min(abs(abs(zz)-refDepth));
  layers_krho = [idx_refDepth 1];    %%% Pressure reference level, level indice k
    % layers_krho = 1 % High-resolution zz (51)=-1.9943e+03 m;
  
  %%% If set true, the GM bolus velocity is added to the calculation
  layers_bolus = false;  
   
  %%% Layers
    for nl=1:layers_maxNum    
      % layers_parm01.addParm(['layers_bounds'],layers_bounds,PARM_REALS); 
      % layers_parm01.addParm(['layers_krho'],layers_krho,PARM_INT); 
      % layers_parm01.addParm(['layers_name'],strtrim(layers_name),PARM_STR); 
      layers_parm01.addParm(['layers_name(' num2str(nl) ')'],strtrim(layers_name(nl,:)),PARM_STR); 
      layers_parm01.addParm(['layers_krho(' num2str(nl) ')'],layers_krho(nl),PARM_INT); 
      layers_parm01.addParm(['layers_bounds(:,' num2str(nl) ')'],layers_bounds(:,nl),PARM_REALS); 
    end
      layers_parm01.addParm('layers_bolus',layers_bolus,PARM_BOOL); 

  
  
  %%z% Create the data.layers file
  write_data_layers(inputpath,LAYERS_PARM,listterm,realfmt);
  
  %%% Create the LAYERS_SIZE.h file
  createLAYERSSIZEh(codepath,length(layers_bounds)-1,layers_maxNum); 
  
 end
  


  

  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
  %%% To store parameter names and values
  diag_parm01 = parmlist;
  diag_parm02 = parmlist;
  DIAG_PARM = {diag_parm01,diag_parm02};
  diag_matlab_parm01 = parmlist;
  DIAG_MATLAB_PARM = {diag_matlab_parm01}; %%% Matlab parameters need to be different
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DIAGNOSTICS PARAMETERS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Stores total number of diagnostic quantities
  ndiags = 0;
       
  diag_parm01.addParm('diag_mnc',false,PARM_BOOL);  
  diag_parm01.addParm('diag_pickup_read',false,PARM_BOOL);  
  diag_parm01.addParm('diag_pickup_write',false,PARM_BOOL);  

  if(run_type~='prod')
        diag_fields_avg = {...   
            %%%%%%%% for spin-up
           'UVEL','WVEL','VVEL','THETA'...
           % 'UVEL','WVEL','VVEL','THETA','UVELSQ','VVELSQ','WVELSQ','THETASQ',...
           % 'Um_Diss','Vm_Diss','Wm_Diss',...
            ... % 'ETAN',...
            ... % 'PHIHYD','PHI_NH',...
            ... % 'DRHODR','SALT',...
        };
  else 
        diag_fields_avg = {...   
            ... %%%%%%%%% for analysis
            'UVEL','VVEL', 'WVEL','THETA','PHIHYD','PHI_NH','DRHODR','ETAN',...
            ... %%% Heat budget
            'TOTTTEND','WTHMASS',...
            'ADVr_TH','ADVx_TH','ADVy_TH','DFxE_TH','DFyE_TH','DFrI_TH','DFrE_TH',...
            'AB_gT','gTinAB',...
            'UVELTH','WVELTH','VVELTH',...
            ... %%% Energy budget
            'momKE','UVELSQ','VVELSQ','WVELSQ',...
            'UV_VEL_Z','WU_VEL','WV_VEL',...
            ...
            'VAHZSMAG','VA4ZSMAG','VAHDSMAG','VA4DSMAG'...
            ...
            'TOTUTEND','TOTVTEND',...
            'Um_Diss','Um_Advec','Um_Ext','Um_AdvZ3','Um_AdvRe','Um_Cori','Um_ImplD','AB_gU',...
            'Vm_Diss','Vm_Advec','Vm_Cori','Vm_Ext','Vm_AdvZ3','Vm_AdvRe','Vm_ImplD','AB_gV',...
            'Wm_Diss','Wm_Advec','WSidDrag','AB_gW',...
            'Um_dPhiX','Vm_dPhiY',... % replace Um_dPHdx and Vm_dPHdy by 'Um_dPhiX' and'Vm_dPhiY' for new code
            'botTauX','botTauY',... %%% Only exist in the new code
            'UVELPHI','VVELPHI',...
            ... %%% KPP diagnostics
            ... %'KPPg_TH','KPPg_SLT','KPPviscA','KPPdiffS','KPPdiffT','KPPghatK','KPPhbl','MXLDEPTH','KPPfrac',...
        };

  end
      
  numdiags_avg = length(diag_fields_avg);  
  diag_freq_avg = 60*t1min;
  % diag_freq_avg = 1*t1day;

  diag_phase_avg = 0;    
      
  for n=1:numdiags_avg    
    ndiags = ndiags + 1;
    diag_parm01.addParm(['fields(1,',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['fileName(',num2str(n),')'],diag_fields_avg{n},PARM_STR);  
    diag_parm01.addParm(['frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_parm01.addParm(['timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL); 
    diag_matlab_parm01.addParm(['diag_fields{1,',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_fileNames{',num2str(n),'}'],diag_fields_avg{n},PARM_STR);  
    diag_matlab_parm01.addParm(['diag_frequency(',num2str(n),')'],diag_freq_avg,PARM_REAL);  
    diag_matlab_parm01.addParm(['diag_timePhase(',num2str(n),')'],diag_phase_avg,PARM_REAL);   
  end
  
if(run_type=='prod')
    diag_fields_inst = {...
        % 'UVEL','VVEL', 'WVEL','THETA',...
            ...%  'UVEL','VVEL', 'WVEL','THETA','PHIHYD','PHI_NH','DRHODR','ETAN',...
            ... % 'RHOAnoma','LaVH1RHO','LaHs1RHO','LaVH2TH','LaHs2TH','LaUH1RHO','LaHw1RHO','LaUH2TH','LaHw2TH',...
            ... % 'ADVr_TH', 'ADVx_TH', 'ADVy_TH','DFrE_TH', 'DFrI_TH', 'DFxE_TH', 'DFyE_TH','TOTTTEND','TRAC01','TRAC02','Tp_gTr01','Tp_gTr02',...
          };
      numdiags_inst = length(diag_fields_inst);  
       % diag_freq_inst = 1*t1day;
      diag_freq_inst = 60*t1min;
      diag_phase_inst = 0;
    
      for n=1:numdiags_inst    
        ndiags = ndiags + 1;
        diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
        diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
        diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
        diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
        diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
        diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
        diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
        diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);     
      end
else
        diag_fields_inst = {...
            %%%%%%%% for spin-up
            % 'UVEL','WVEL','THETA'...
            ... % 'DRHODR','SALT','PHIHYD','PHI_NH','VVEL','ETAN'...
          };
      numdiags_inst = length(diag_fields_inst);  
      diag_freq_inst = 60*t1min;
      % diag_freq_inst = 10;
      diag_phase_inst = 0;
    
      for n=1:numdiags_inst    
        ndiags = ndiags + 1;
        diag_parm01.addParm(['fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
        diag_parm01.addParm(['fileName(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
        diag_parm01.addParm(['frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
        diag_parm01.addParm(['timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL); 
        diag_matlab_parm01.addParm(['diag_fields(1,',num2str(ndiags),')'],diag_fields_inst{n},PARM_STR);  
        diag_matlab_parm01.addParm(['diag_fileNames(',num2str(ndiags),')'],[diag_fields_inst{n},'_inst'],PARM_STR);  
        diag_matlab_parm01.addParm(['diag_frequency(',num2str(ndiags),')'],-diag_freq_inst,PARM_REAL);  
        diag_matlab_parm01.addParm(['diag_timePhase(',num2str(ndiags),')'],diag_phase_inst,PARM_REAL);     
      end
end
  

  %%% Create the data.diagnostics file
  write_data_diagnostics(inputpath,DIAG_PARM,listterm,realfmt);
  
  %%% Create the DIAGNOSTICS_SIZE.h file
  if(useLAYERS)
    createDIAGSIZEh(codepath,ndiags,max(Nr,length(layers_bounds)-1));
  else
    createDIAGSIZEh(codepath,ndiags,Nr);
  end
 


  if(useOBCS)
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
  %%%%% OBCS %%%%%
  %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%
    
  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OBCS SET-UP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% To store parameter names and values
  %%% Add 2-element cell arrays to this cell array in the form 
  %%%  OBCS_PARM{1} = addParameter(OBCS_PARM{1},'paramName',paramValue,parmType);
  %%% to specify additional parameters. The parameter type parmType must
  %%% take one of the integer values above.
    %%% To store parameter names and values
  obcs_parm01 = parmlist;
  obcs_parm02 = parmlist;
  obcs_parm03 = parmlist;
  obcs_parm04 = parmlist;
  obcs_parm05 = parmlist;
  OBCS_PARM = {obcs_parm01,obcs_parm02,obcs_parm03,obcs_parm05};  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% DEFINE OPEN BOUNDARY TYPES (OBCS_PARM01) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Set boundary points that are open   
  if (useobcsNorth)
      OB_Jnorth= Ny*ones(1,Nx);
      obcs_parm01.addParm('OB_Jnorth',OB_Jnorth,PARM_INTS); 
  end
  if(useobcsSouth)
      OB_Jsouth = ones(1,Nx);
      obcs_parm01.addParm('OB_Jsouth',OB_Jsouth,PARM_INTS); 
  end
  
  % tidalPeriod = 43200;
  tidalPeriod= [43200,43200,43200,43200,43200,43200,43200,43200,43200,43200];
  % tidalPeriod= [86400,86400,86400,86400,86400,86400,86400,86400,86400,86400];
%   tidalPeriod=[44714.16,43200.,45569.88,43081.92,86164.2,92949.48,86637.24,96726.24,1180295.64,2380716];
  obcs_parm01.addParm('useOBCStides',useOBCStides,PARM_BOOL);
  
  if (useOBCStides)
          obcs_parm01.addParm('tidalPeriod',tidalPeriod,PARM_INTS);    
     if (useobcsNorth)
          OBNamFile= 'OBNam.obcs';
          OBNphFile= 'OBNph.obcs';
          obcs_parm01.addParm('OBNamFile',OBNamFile,PARM_STR);  
          obcs_parm01.addParm('OBNphFile',OBNphFile,PARM_STR); 
     end
     if (useobcsSouth)
          OBSamFile= 'OBSam.obcs';
          OBSphFile= 'OBSph.obcs';
          obcs_parm01.addParm('OBSamFile',OBSamFile,PARM_STR);  
          obcs_parm01.addParm('OBSphFile',OBSphFile,PARM_STR); 
     end

     % create tidal input files
         tidalComponents=10;
         OBns ={};
         if (useobcsNorth)
             OBns = {'N'};
         end
         if (useobcsSouth)
             OBns = [OBns, {'S'}];
         end
     
        for ob = OBns
            OBlength=Ny;
            if any(strcmp(ob,{'N','S'}))
                OBlength=Nx;
            end
            for fld={'am','ph'}
                fnm=['OB' ob{1} fld{1} '.obcs'];
                % tmp=randn(OBlength,tidalComponents)/1000;
                tmp=zeros(OBlength,tidalComponents);

                 % Phase = 0;
                 Phase = 3 * 3600;
                 Deep_lead_Shallow = 0;
                 % Deep_lead_Shallow = 1*3600; %%% From Kurt: The tides in the deep ocean lead the tide in the shallow ocean for about 30 degrees, i.e., about 1 hour

                % specify (0.1 m/s, 2 hr) for North boundary tidal component 1
                if strcmp(ob,'N')
                    if strcmp(fld,'am')
                        tmp(:,1) = tmp(:,1) + Atide;
                    else
                        if(varyingtidalphase)
                            varyingphase = 0.5*t1hour;
                            for iPH = 1:Nx
                                tmp(iPH,1)=tmp(iPH,1)+(Phase-varyingphase/2)+varyingphase/Nx*iPH;
                            end
                        else
                            tmp(:,1) = tmp(:,1) + Phase;
                        end
                    end
                end
                % specify (0.1 m/s, 2 hr) for South boundary tidal component 1
                if strcmp(ob,'S')
                    if strcmp(fld,'am')
                        Atide_south = Atide*h(end)/h(1);
                        tmp(:,1) = tmp(:,1) + Atide_south;
                    else
                        if(varyingtidalphase)
                            varyingphase = 0.5*t1hour;
                            for iPH = 1:Nx
                                tmp(iPH,1)=tmp(iPH,1)+(Phase-varyingphase/2)+varyingphase/Nx*iPH;
                            end
                        else
                            tmp(:,1) = tmp(:,1) + Phase - Deep_lead_Shallow;
                        end
                    end
                end
                
                writeDataset(tmp,fullfile(inputpath,fnm),ieee,prec);
            end
        end

  else 
      Atide = 0;
  end
  
 
  
  %%% Enforces mass conservation across the northern boundary by adding a
  %%% barotropic inflow/outflow  
  useOBCSbalance = false;  
  obcs_parm01.addParm('useOBCSbalance',useOBCSbalance,PARM_BOOL);

  if(useOBCSbalance)
        OBCS_balanceFacN = 1; %%% A value -1 balances an individual boundary
        OBCS_balanceFacS = 1;
        obcs_parm01.addParm('OBCS_balanceFacN',OBCS_balanceFacN,PARM_REAL); 
        obcs_parm01.addParm('OBCS_balanceFacS',OBCS_balanceFacS,PARM_REAL);  
  end
  %%% Enables/disables sponge layers   
  useOBCSsponge = true;
  obcs_parm01.addParm('useOBCSsponge',useOBCSsponge,PARM_BOOL);

  %%% Set boundary velocities and temperatures to be consistent with the
  %%% streamfunction psi = alpha*x*y, u=-dpsi/dy, v=dpsi/dx
  useOBCSprescribe = true;  
  
  OBNt = ones(Nx,1)*tNorth;
  OBNs = ones(Nx,1)*sNorth;
  OBNv = ones(Nx,1)*vNorth;
  OBSt = ones(Nx,1)*tSouth;
  OBSs = ones(Nx,1)*sSouth;
  OBSv = ones(Nx,1)*vSouth;

  %%% Write boundary variables to files  
  %%% Set OBCS prescription parameters
  obcs_parm01.addParm('useOBCSprescribe',useOBCSprescribe,PARM_BOOL);
  if (useobcsNorth)
      writeDataset(OBNt,fullfile(inputpath,'OBNtFile.bin'),ieee,prec);
      writeDataset(OBNs,fullfile(inputpath,'OBNsFile.bin'),ieee,prec);
      writeDataset(OBNv,fullfile(inputpath,'OBNvFile.bin'),ieee,prec);
      obcs_parm01.addParm('OBNtFile','OBNtFile.bin',PARM_STR);
      obcs_parm01.addParm('OBNsFile','OBNsFile.bin',PARM_STR);
      obcs_parm01.addParm('OBNvFile','OBNvFile.bin',PARM_STR);
  end
  if (useobcsSouth)
      writeDataset(OBSt,fullfile(inputpath,'OBStFile.bin'),ieee,prec);
      writeDataset(OBSs,fullfile(inputpath,'OBSsFile.bin'),ieee,prec);
      writeDataset(OBSv,fullfile(inputpath,'OBSvFile.bin'),ieee,prec);
      obcs_parm01.addParm('OBStFile','OBStFile.bin',PARM_STR);
      obcs_parm01.addParm('OBSsFile','OBSsFile.bin',PARM_STR);
      obcs_parm01.addParm('OBSvFile','OBSvFile.bin',PARM_STR);
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ORLANSKI OPTIONS (OBCS_PARM02) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Enables/disables Orlanski radiation conditions at the boundaries -
  %%% allows waves to propagate out through the boundary with minimal
  %%% reflection  
  obcs_parm01.addParm('useOrlanskiNorth',useOrlanskiNorth,PARM_BOOL);
  obcs_parm01.addParm('useOrlanskiSouth',useOrlanskiSouth,PARM_BOOL);
  %%% Velocity averaging time scale - must be larger than deltaT.
  %%% The Orlanski radiation condition computes the characteristic velocity
  %%% at the boundary by averaging the spatial derivative normal to the 
  %%% boundary divided by the time step over this period.
  %%% At the moment we're using the magic engineering factor of 3.
  cvelTimeScale = 3*deltaT; % Averaging period for phase speed (s)
  %%% Max dimensionless CFL for Adams-Basthforth 2nd-order method
  CMAX = 0.45; 
  
  if(useOrlanskiNorth || useOrlanskiSouth)
  obcs_parm02.addParm('cvelTimeScale',cvelTimeScale,PARM_REAL);
  obcs_parm02.addParm('CMAX',CMAX,PARM_REAL);
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% SPONGE LAYER OPTIONS (OBCS_PARM03) %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   obcs_parm03.addParm('spongeThickness',spongeThickness,PARM_INT);

    Vrelaxobcsinner = 12*t1hour;
    Vrelaxobcsbound = 1*t1hour;
    % Vrelaxobcsinner = 3*t1hour;
    % Vrelaxobcsbound = 1*t1min;
%   %% Relaxation time at meridional boundaries set to time for inflow to
%   %% cross the sponge layer
%   Vrelaxobcsbound = spongeThicknessDim/(abs(alpha)*Ly/2);
  
  obcs_parm03.addParm('Vrelaxobcsinner',Vrelaxobcsinner,PARM_REAL);
  obcs_parm03.addParm('Vrelaxobcsbound',Vrelaxobcsbound,PARM_REAL);

  if(useobcsEast||useobcsWest)
        Urelaxobcsinner = 864000;
        Urelaxobcsbound = 43200;
        obcs_parm03.addParm('Urelaxobcsinner',Urelaxobcsinner,PARM_REAL);
        obcs_parm03.addParm('Urelaxobcsbound',Urelaxobcsbound,PARM_REAL);  
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE THE 'data.obcs' FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  %%% Creates the 'data.obcs' file
  write_data_obcs(inputpath,OBCS_PARM,listterm,realfmt);

  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%% PACKAGES %%%%%
  %%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%
  
  packages = parmlist;
  PACKAGE_PARM = {packages};  
  
  packages.addParm('useDiagnostics',true,PARM_BOOL);    
  packages.addParm('useKPP',false,PARM_BOOL);
  packages.addParm('useRBCS',useRBCS,PARM_BOOL);        
  packages.addParm('useEXF',useEXF,PARM_BOOL);        
  packages.addParm('useCAL',useEXF,PARM_BOOL); 
  packages.addParm('useOBCS',useOBCS,PARM_BOOL);  
  packages.addParm('useLAYERS',useLAYERS,PARM_BOOL);  

  %%% Create the data.pkg file
  write_data_pkg(inputpath,PACKAGE_PARM,listterm,realfmt);
  
  
 
   



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WRITE PARAMETERS TO A MATLAB FILE %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%% Creates a matlab file defining all input parameters
  ALL_PARMS =[PARM PACKAGE_PARM DIAG_MATLAB_PARM];
  if (useRBCS)
      ALL_PARMS = [ALL_PARMS RBCS_PARM];
  end
  if (useEXF)
      ALL_PARMS = [ALL_PARMS EXF_PARM];
  end
  if (useLAYERS)
    ALL_PARMS = [ALL_PARMS LAYERS_PARM];
  end  
  %%% Creates a matlab file defining all input parameters
  write_matlab_params(inputpath,ALL_PARMS,realfmt);
  
  
end



    %%% Specifies shape of coastal walls. Must satisfy f=1 at x=0 and f=0 at
    %%% x=1.
    %%%
    function f = coastShape (x)
     
      f = 0.5.*(1+cos(pi*x));
    %   f = exp(-x);
    %   f = 1-x;
      
    end

