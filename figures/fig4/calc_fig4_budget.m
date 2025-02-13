%%% Calculate and store the buoyancy budget for fig3
load('fig4/fig4_km.mat')

tidx = nt_percycle*5+1:nt_percycle*10;
tt = tt(tidx);


%--- Compute the vorticity budget
re_zeta = real(zeta);
dzetadt = [0 (re_zeta(3:end)-re_zeta(1:end-2))/dt/2 0];
dzetadt = dzetadt(tidx);

bxcs = real(1i*kx.*buoy(tidx)*cs);
bzss = real(-1i*mz_t(tidx).*buoy(tidx)*ss);

%--- Normalization
bxcs = bxcs/max(abs(dzetadt));
bzss = bzss/max(abs(dzetadt));
dzetadt = dzetadt/max(abs(dzetadt));


%--- make a plot
figure(2)
clf;set(gcf,'Color','w')
plot(tt/43200,bxcs,'LineWidth',2)
hold on;
plot(tt/43200,bzss,'LineWidth',2)
plot(tt/43200,dzetadt,'--','LineWidth',2)
plot(tt/43200,dzetadt-bxcs-bzss,':','LineWidth',2)
xlim([5 10])
set(gca,'Fontsize',16)
grid on;grid minor;




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
