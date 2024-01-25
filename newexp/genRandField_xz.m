%%%
%%% genRandField_xz.m
%%%
%%% Generates a random two-dimensional field. Useful for creating random 
%%% initial conditions or random topography. The output field 'psi' will
%%% have characteristic spectral wavelength 'lambda' with an exponential 
%%% width 'W' in spectral space, and and rms amplitude 'psirms'.
%%% Nx and Nr define the grid size, and Lx and H define the 
%%% zonal and vertial domain lengths. If 'W' is an empty vector then a
%%% default value of 1/8 * the wavenumber corresponding to lambda will be
%%% used.
%%% 
function F = genRandField (lambda,W,Frms,Nx,Nr,Lx,H) 
 
  %%% Spectral grids  
  k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
  K_xk = 2*pi.*(k)./Lx;
  l = [0:1:Nr/2-1,-Nr/2:1:-1]; %%% Meridional wavenumber
  K_zl = 2*pi.*(l)./H;
  [K_zkl,K_xkl] = meshgrid(K_zl, K_xk);   
  
  %%% Most energetic wavenumber
  K_0 = 2*pi/lambda; 
  
  %%% Exponential width of energy band in wavenumber space
  if (isempty(W))
    W = K_0/8; 
  end

  %%% Amplitude is exponential about K0, and phases are random. N.B. here
  %%% we only define the amplitude up to a constant - below we constrain it.  
  K = sqrt(K_xkl.^2 + K_zkl.^2);
  theta = 2 .* pi .* rand(Nx,Nr);
  psi_fft = K.^(-1).*exp(-((K-K_0)/W).^2) .* exp(1i*theta);

  %%% Avoids infinite mode-0 amplitude 
  psi_fft(1,1) = 0;

  %%% Transform back to real space
  F = real(ifft2(psi_fft));

  %%% Normalize so that RMS is Frms
  F = F * Frms./sqrt(sum(sum(F.^2))/(Nx*Nr));  

end