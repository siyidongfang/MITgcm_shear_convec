
figure(1)
pcolor(re_d2bdz2)
shading flat;colorbar;colormap(redblue);clim([-500 500]/10)

figure(2)
pcolor(re_dbdz)
shading flat;colorbar;colormap(redblue);clim([-5 5])


aaa = re_dbdz(:,1);
bbb = re_dbdz(:,Nr+1);
figure(3)
clf;
plot(aaa);
hold on;
plot(bbb,'--')


aaa = re_d2bdz2(:,1);
bbb = re_d2bdz2(:,Nr);
figure(4)
clf;
plot(aaa);
hold on;
plot(bbb,'--')


aaa = re_d2zetadz2(:,1);
bbb = re_d2zetadz2(:,Nr+1);
figure(5)
clf;
plot(aaa);
hold on;
plot(bbb,'--')


figure(6)
pcolor(re_d2zetadz2)
shading flat;colorbar;colormap(redblue);
clim([-5 5]/1e16)
