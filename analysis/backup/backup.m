%%% Calculate imposed large-scale velocity
fid = fopen(fullfile(exppath,'/input','uVelInitFile.bin'),'r','b');
Az = zeros(Nx,Ny,Nr);
for k=1:Nr
  Az(:,:,k) = fread(fid,[Nx Ny],'real*8');
end
fclose(fid);
Az = squeeze(Az(1,:));
omega = 2*pi/43200;
DT = 3600;
topo=0;

u12 = zeros(12,Nr);
v12 = zeros(12,Nr);
for o=1:12
    t1 = (o-1)*3600;
    t2 = o*3600;
    u12(o,:) = Az/omega/DT*(sin(omega*t2)-sin(omega*t1));
    v12(o,:) = Az*f0*cosd(topo)/(omega^2)/DT*(cos(omega*t2)-cos(omega*t1));
end
aa = [u12(:,end);u12(:,end)];
bb = [v12(:,end);v12(:,end)];
figure(1);plot(aa);hold on;plot(bb);hold off;
