% load('grow_K1.mat')
load('grow_M2.mat')

Ri_flat = Ri(:);
[Ri_sort,I] = sort(Ri_flat);

grow_flat = grow_floquet(:);
grow_sort = grow_flat(I);

isConvec_flat = isConvec(:);
isConvec_sort = isConvec_flat(I);

grow_sort(isConvec_sort==1)=NaN;
Ri_sort(isConvec_sort==1)=NaN;

grow_sort(grow_sort<=0)=NaN;


Ri_sort = Ri_sort(~isnan(grow_sort));
grow_sort = grow_sort(~isnan(grow_sort));


fontsize = 18;
figure(1);set(gcf,'Color','w')
scatter(1./Ri_sort,grow_sort)
xlim([0 10])
ylim([0 0.75])
grid on;box on;
set(gca,'Fontsize',fontsize)

