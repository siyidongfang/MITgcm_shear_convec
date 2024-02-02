
clear;
Ny = 1;
Nx = 150;
Nr = 318;
Lx = 3000;
H = 1500;
tNoise = 1e-6;
noise_length = 10;
noise_amp = 1;
W = [];
% W = 0.01;


%%% Varied dz with depth  %  -- from Xiaozhou
Hmax = 1500;
Hsurface = 1000;
Ntop = 150;
dz_const = 3;
dz = dz_const.*ones(1,Nr);
dz(Nr-Ntop + 1:Nr) = dz(Nr - Ntop) * 1.009.^(1:Ntop);
sum_dz_sponge = sum(dz(Nr-Ntop + 1:Nr));
dz(Nr-Ntop + 1:Nr) = dz(Nr-Ntop + 1:Nr).*Hsurface/sum_dz_sponge;
dz = flipud(dz')';
zz = -cumsum((dz+[0 dz(1:end-1)])/2);

dx = Lx/Nx*ones(1,Nx);
xx = cumsum((dx + [0 dx(1:end-1)])/2);

Nx_noise = Lx;
Nr_noise = Hmax;

Fnoise = genRandField_xz(noise_length,W,noise_amp,Nx_noise,Nr_noise,Lx,Hmax);
Fnoise = tNoise*Fnoise/max(max(abs(Fnoise)));

[zzz1,xxx1] = meshgrid(-1*(1:Hmax),1:Lx);
[zzz2,xxx2] = meshgrid(zz,xx);

Fnoise_interp = interp2(zzz1,xxx1,Fnoise,zzz2,xxx2);

% midpointX = round(Nx/2);
% Fnoise2 = [Fnoise(midpointX+1:Nx,:);Fnoise(1:midpointX,:)];

figure(1)
pcolor(Fnoise');shading flat;colorbar;colormap(redblue);
clim([-1 1]*tNoise)



figure(2)
clf;set(gcf,'Color','w')
pcolor(xx,-zz,Fnoise_interp');shading flat;colorbar;colormap(redblue);
clim([-1 1]*tNoise);axis ij;
fontsize = 17;
xlabel('Distance, y (km)');
ylabel('Depth (m)');
title({'Initial temperature (^oC)',['Characteristic spectral wavelength \lambda = ' num2str(noise_length) 'm']});
set(gca,'fontsize',fontsize);
PLOT = gcf;
PLOT.Position =[80 341 606 343];
grid on;grid minor;

% %%
% hydroTh = ones(Nx,Ny,Nr);
% for i=1:Nx
%     for j=1:Ny
%         for k=1:Nr
%             hydroTh(i,j,k)=Fnoise(i,k);
%         end
%     end
% end

        


