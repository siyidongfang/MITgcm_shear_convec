%%%
%%% genRandField_y.m
%%%
%%% Generates a random 1-dimensional field. Useful for creating random 
%%% initial conditions or random topography. The output field 'psi' will
%%% have characteristic spectral wavelength 'lambda' with an exponential 
%%% width 'W' in spectral space, and and rms amplitude 'psirms'.
%%% Ny define the grid size, and Ly define the 
%%% meridional domain lengths. If 'W' is an empty vector then a
%%% default value of 1/8 * the wavenumber corresponding to lambda will be
%%% used.
%%% 
function F = genRandField_y (lambda,W,Frms,Ny,Ly) 
 
  %%% Spectral grids  
  l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
  K_yl = 2*pi.*(l)./Ly;  
  K_ykl = K_yl;

  %%% Most energetic wavenumber
  K_0 = 2*pi/lambda; 
  
  %%% Exponential width of energy band in wavenumber space
  if (isempty(W))
    W = K_0/8; 
  end

  %%% Amplitude is exponential about K0, and phases are random. N.B. here
  %%% we only define the amplitude up to a constant - below we constrain it.  
  K = K_ykl;
  % theta = 2 .* pi .* rand(1,Ny-1);
  theta = 2 .* pi .* rand(1,Ny);
  psi_fft = K.^(-1).*exp(-((K-K_0)/W).^2) .* exp(1i*theta);

  %%% Avoids infinite mode-0 amplitude 
  psi_fft(1,1) = 0;

  %%% Transform back to real space
  F = real(ifft(psi_fft));

  %%% Normalize so that RMS is Frms
  F = F * Frms./sqrt(sum(F.^2)/Ny);  

end