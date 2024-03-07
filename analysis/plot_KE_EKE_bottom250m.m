%%%
%%% plot_KE_EKE_T_S_series.m
%%%
%%% Plots the total kinetic energy output from MITgcm simulations.
%%%
%%% NOTE: Doesn't account for u/v gridpoint locations, and doesn't handle
%%% partial cells.
%%%

    
clear; close all;

ne=4;
load_all;

%%% Add path
addpath functions/
addpath colormaps;
addpath colormaps/cmocean/;

figdir = [exppath '/img/'];

% nDumps = 120

ntime = zeros(1,nDumps);
EKEtot = zeros(1,nDumps);
KEtot = zeros(1,nDumps);
KElen = 0;

nseries = 1:nDumps;

for n= nseries
 
  ntime(n) =  dumpIters(n)*deltaT/86400;  
  
  uvel = rdmdsWrapper(fullfile(exppath,'/results/UVEL'),dumpIters(n));      
  vvel = rdmdsWrapper(fullfile(exppath,'/results/VVEL'),dumpIters(n));      
  wvel = rdmdsWrapper(fullfile(exppath,'/results/WVEL'),dumpIters(n));      
  uvelsq = rdmdsWrapper(fullfile(exppath,'/results/UVELSQ'),dumpIters(n));      
  vvelsq = rdmdsWrapper(fullfile(exppath,'/results/VVELSQ'),dumpIters(n));      
  wvelsq = rdmdsWrapper(fullfile(exppath,'/results/WVELSQ'),dumpIters(n));      
  
  Nbot = round(250/delR(end));
  uvel = uvel(:,:,Nr-Nbot:Nr);
  vvel = vvel(:,:,Nr-Nbot:Nr);
  wvel = wvel(:,:,Nr-Nbot:Nr);
  uvelsq = uvelsq(:,:,Nr-Nbot:Nr);
  vvelsq = vvelsq(:,:,Nr-Nbot:Nr);
  wvelsq = wvelsq(:,:,Nr-Nbot:Nr);

  if (isempty(uvelsq) || isempty(vvelsq) || isempty(wvelsq))
    break;
  end

  KE = 0.5*(uvelsq + vvelsq + wvelsq);
  EKE = 0.5*(uvelsq + vvelsq + wvelsq - uvel.^2 - vvel.^2 - wvel.^2);
  KEtot(n) = 0;
  EKEtot(n) = 0;
  

  for i=1:Nx
    for j=1:Ny
      for k=1:Nbot+1
        KEtot(n) = KEtot(n) + KE(i,j,k)*delX(i)*delY(j)*delR(k);
        EKEtot(n) = EKEtot(n) + EKE(i,j,k)*delX(i)*delY(j)*delR(k);
      end
    end
  end
  KElen = KElen + 1;
end
  

DX_xyz = repmat(reshape(delX,[Nx 1 1]),[1 Ny Nr]);
DY_xyz = repmat(reshape(delY,[1 Ny 1]),[Nx 1 Nr]);
DZ_xyz = repmat(reshape(delR,[1 1 Nr]),[Nx Ny 1]);
Tseries = zeros(1,nDumps);
% Sseries = zeros(1,nDumps);


% for n=1:KElen
for n=nseries
  theta = rdmdsWrapper(fullfile(exppath,'/results/THETA'),dumpIters(n));  
  % salt  = rdmdsWrapper(fullfile(exppath,'/results/SALT') ,dumpIters(n));  
  Vol = sum(hFacC(:,:,Nr-Nbot:Nr).*DX_xyz(:,:,Nr-Nbot:Nr).*DY_xyz(:,:,Nr-Nbot:Nr).*DZ_xyz(:,:,Nr-Nbot:Nr),'all');

  Tseries(n) = sum(theta(:,:,Nr-Nbot:Nr).*hFacC(:,:,Nr-Nbot:Nr).*DX_xyz(:,:,Nr-Nbot:Nr).*DY_xyz(:,:,Nr-Nbot:Nr).*DZ_xyz(:,:,Nr-Nbot:Nr),'all')/Vol;
  % Sseries(n) = sum(salt.*hFacC.*DX_xyz.*DY_xyz.*DZ_xyz,'all')/Vol;  
end



%%
%%% Make plots!
fontsize = 16;

figure(1);
clf;
set(gcf,'color','w')
set(gcf,'Position',[61 136 1337 814])
subplot(3,1,1)
plot(ntime(nseries)*2,KEtot(nseries),'LineWidth',2);
axis tight;
xlabel('Tidal cycles','interpreter','latex');
ylabel('Total KE (m^2s^-^2)');
set(gca,'FontSize',fontsize);
grid on;grid minor;
subplot(3,1,2)
plot(ntime(nseries)*2,EKEtot(nseries),'LineWidth',2);
axis tight;
xlabel('Tidal cycles','interpreter','latex');
ylabel('TKE (m^2s^-^2)');
set(gca,'FontSize',fontsize);
grid on;grid minor;
subplot(3,1,3)
plot(ntime(nseries)*2,Tseries(nseries),'LineWidth',2);
axis tight;
xlabel('Tidal cycles','interpreter','latex');
ylabel('\theta (degC)');
set(gca,'FontSize',fontsize);
grid on;grid minor;
% subplot(2,2,4)
% plot(ntime(nseries)*2,Sseries(nseries),'LineWidth',2);
% axis tight;
% xlabel('Time (days)');
% ylabel('S (psu)');
% grid on;grid minor;
% set(gca,'FontSize',fontsize);
% print('-dpng','-r150',[figdir 'series_KE_EKE_T_S.png']);
