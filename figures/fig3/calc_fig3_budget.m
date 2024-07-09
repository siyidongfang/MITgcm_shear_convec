%%% Calculate and store the buoyancy budget for fig3
load('fig3/fig3_km.mat')

tidx = nt_percycle*5+1:nt_percycle*10;

%--- Calculate the buoyancy budget
uB0x = -re_uuu(tidx)*N^2*ss;
wB0z = -re_www(tidx)*N^2*cs;
wBz  =  re_www(tidx)*shear/omega*N^2*ss.*st(tidx);
% diffusion = 0*
dbdt = [0 (re_buoy(3:end)-re_buoy(1:end-2))/dt/2 0];
dbdt = dbdt(tidx);

%--- Normalization
uB0x = uB0x/max(abs(dbdt));
wB0z = wB0z/max(abs(dbdt));
wBz = wBz/max(abs(dbdt));
dbdt = dbdt/max(abs(dbdt));
tt = tt(tidx);

%--- make a plot
figure(1)
clf
plot(tt/43200,uB0x)
hold on;
plot(tt/43200,wB0z)
plot(tt/43200,wBz)
plot(tt/43200,dbdt,'--')
plot(tt/43200,dbdt-uB0x-wB0z-wBz,':')
xlim([5 10])

%--- save the data