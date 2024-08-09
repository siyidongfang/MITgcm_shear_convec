clear all;close all;

addpath ../instability_eig/

load('grow_M2_calc.mat')
mycolormap = WhiteBlueGreenYellowRed(8);
fontsize = 18;markersize = 36;

figure(1);clf;
scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1200 800],'Color','w');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1450 750],'Color','w');
Nidx = 1:Nn;

for nt =Ntopo:-1:1
    topo = topo_all(nt);

    Ri_flat = Ri(nt,Nidx,:); Ri_flat = Ri_flat(:);
    grow_flat = grow_floquet(nt,Nidx,:); grow_flat = grow_flat(:);
    Bu_flat = Bu(nt,Nidx,:); Bu_flat = Bu_flat(:);
    strat_flat = strat(nt,Nidx,:); strat_flat = strat_flat(:);
    isConvec_flat = isConvec(nt,Nidx,:); isConvec_flat = isConvec_flat(:);
    
    [Ri_sort,I] = sort(Ri_flat);
    grow_sort = grow_flat(I);
    Bu_sort = Bu_flat(I);
    strat_sort = strat_flat(I);
    isConvec_sort = isConvec_flat(I);
    
    grow_sort(isConvec_sort==1)=NaN;
    Ri_sort(isConvec_sort==1)=NaN;

    Ri_sort = Ri_sort(~isnan(grow_sort));
    grow_sort = grow_sort(~isnan(grow_sort));
    Bu_sort = Bu_sort(~isnan(grow_sort));
    strat_sort = strat_sort(~isnan(grow_sort)); 

    subplot(2,2,1);hold on;
    scatter(Ri_sort,grow_sort,markersize,topo.*ones(1,length(grow_sort)),'filled')
end



for nn =Nidx
    N = N_all(nn);

    Ri_flat = Ri(:,nn,:); Ri_flat = Ri_flat(:);
    grow_flat = grow_floquet(:,nn,:); grow_flat = grow_flat(:);
    Bu_flat = Bu(:,nn,:); Bu_flat = Bu_flat(:);
    isConvec_flat = isConvec(:,nn,:); isConvec_flat = isConvec_flat(:);
    
    [Ri_sort,I] = sort(Ri_flat);
    grow_sort = grow_flat(I);
    Bu_sort = Bu_flat(I);
    isConvec_sort = isConvec_flat(I);
    
    grow_sort(isConvec_sort==1)=NaN;
    Ri_sort(isConvec_sort==1)=NaN;

    Ri_sort = Ri_sort(~isnan(grow_sort));
    grow_sort = grow_sort(~isnan(grow_sort));
    Bu_sort = Bu_sort(~isnan(grow_sort));
    
    subplot(2,2,2);hold on;
    scatter(Ri_sort,grow_sort,markersize,N.*ones(1,length(grow_sort)),'filled')
end


% Colored by Shear
velocityshear_flat = velocityshear(:);
[shear_sort,I] = sort(velocityshear_flat);

Ri_flat = Ri(:);
grow_flat = grow_floquet(:);
Bu_flat = Bu(:);
isConvec_flat = isConvec(:);

Ri_sort = Ri_flat(I);
grow_sort = grow_flat(I);
Bu_sort = Bu_flat(I);
isConvec_sort = isConvec_flat(I);

grow_sort(isConvec_sort==1)=NaN;
Ri_sort = Ri_sort(~isnan(grow_sort));
Bu_sort = Bu_sort(~isnan(grow_sort));
shear_sort = shear_sort(~isnan(grow_sort));
grow_sort = grow_sort(~isnan(grow_sort));

subplot(2,2,3);
scatter(Ri_sort,grow_sort,markersize,log10(shear_sort),'filled')


% Colored by Slope Burger Number
subplot(2,2,4);
scatter(Ri_sort,grow_sort,markersize,(Bu_sort),'filled')


%%


load('grow_K1_calc.mat')

Nidx = 1:Nn;

for nt =Ntopo:-1:1
% for nt =2
    topo = topo_all(nt);

    Ri_flat = Ri(nt,Nidx,:); Ri_flat = Ri_flat(:);
    grow_flat = grow_floquet(nt,Nidx,:); grow_flat = grow_flat(:);
    Bu_flat = Bu(nt,Nidx,:); Bu_flat = Bu_flat(:);
    strat_flat = strat(nt,Nidx,:); strat_flat = strat_flat(:);
    isConvec_flat = isConvec(nt,Nidx,:); isConvec_flat = isConvec_flat(:);
    
    [Ri_sort,I] = sort(Ri_flat);
    grow_sort = grow_flat(I);
    Bu_sort = Bu_flat(I);
    strat_sort = strat_flat(I);
    isConvec_sort = isConvec_flat(I);
    
    grow_sort(isConvec_sort==1)=NaN;
    Ri_sort(isConvec_sort==1)=NaN;
    % grow_sort(grow_sort<=0)=NaN;

    % grow_sort(Ri_sort>100)=NaN;

    Ri_sort = Ri_sort(~isnan(grow_sort));
    grow_sort = grow_sort(~isnan(grow_sort));
    Bu_sort = Bu_sort(~isnan(grow_sort));
    strat_sort = strat_sort(~isnan(grow_sort)); 
    
    % subplot(2,2,1);hold on;
    % if(topo<=10)
    %     scatter(1./Ri_sort,grow_sort,markersize,topo.*ones(1,length(grow_sort)),'filled')
    % else
    %     scatter(1./Ri_sort,grow_sort,markersize,topo.*ones(1,length(grow_sort)))    
    % end
    % xlim([0 10])
    % ylim([0 0.85])
    % grid on;box on;
    % set(gca,'Fontsize',fontsize)
    % colormap(mycolormap);colorbar;
    % clim([0 20]);
    % xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
    % ylabel('(hour$^{-1}$)','interpreter','latex');

    subplot(2,2,1);hold on;
    scatter(Ri_sort,grow_sort,markersize,topo.*ones(1,length(grow_sort)))
    set(gca,'xscale','log')
    xlim([0 120])
    ylim([0 0.85])
    grid on;box on;
    set(gca,'Fontsize',fontsize)
    colormap(mycolormap);colorbar;
    clim([0 20]);
    % xlabel('Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
    ylabel('(hour$^{-1}$)','interpreter','latex');

    title('Growth Rate vs. ${R_i}_\mathrm{min}$ Colored by Slope, $\theta$ ($^\circ$)','interpreter','latex');

end



for nn =Nidx
    N = N_all(nn);

    Ri_flat = Ri(:,nn,:); Ri_flat = Ri_flat(:);
    grow_flat = grow_floquet(:,nn,:); grow_flat = grow_flat(:);
    Bu_flat = Bu(:,nn,:); Bu_flat = Bu_flat(:);
    isConvec_flat = isConvec(:,nn,:); isConvec_flat = isConvec_flat(:);
    
    [Ri_sort,I] = sort(Ri_flat);
    grow_sort = grow_flat(I);
    Bu_sort = Bu_flat(I);
    isConvec_sort = isConvec_flat(I);
    
    grow_sort(isConvec_sort==1)=NaN;
    Ri_sort(isConvec_sort==1)=NaN;
    % grow_sort(grow_sort<=0)=NaN;

    % grow_sort(Ri_sort>100)=NaN;

    Ri_sort = Ri_sort(~isnan(grow_sort));
    grow_sort = grow_sort(~isnan(grow_sort));
    Bu_sort = Bu_sort(~isnan(grow_sort));
    
    subplot(2,2,2);hold on;
    scatter(Ri_sort,grow_sort,markersize,N.*ones(1,length(grow_sort)))
    set(gca,'xscale','log')
    xlim([0 120])
    ylim([0 0.85])
    grid on;box on;
    set(gca,'Fontsize',fontsize)
    colormap(mycolormap);colorbar;
    % xlabel('Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
    ylabel('(hour$^{-1}$)','interpreter','latex');
    title('Growth Rate vs. ${R_i}_\mathrm{min}$ Colored by $\tilde N$ (s$^{-1}$)','interpreter','latex');

end


% Colored by Shear
velocityshear_flat = velocityshear(:);
[shear_sort,I] = sort(velocityshear_flat);

Ri_flat = Ri(:);
grow_flat = grow_floquet(:);
Bu_flat = Bu(:);
isConvec_flat = isConvec(:);

Ri_sort = Ri_flat(I);
grow_sort = grow_flat(I);
Bu_sort = Bu_flat(I);
isConvec_sort = isConvec_flat(I);

grow_sort(isConvec_sort==1)=NaN;
Ri_sort = Ri_sort(~isnan(grow_sort));
Bu_sort = Bu_sort(~isnan(grow_sort));
shear_sort = shear_sort(~isnan(grow_sort));
grow_sort = grow_sort(~isnan(grow_sort));


subplot(2,2,3);hold on;
scatter(Ri_sort,grow_sort,markersize,log10(shear_sort))
set(gca,'xscale','log')
xlim([0 120])
ylim([0 0.85])
grid on;box on;
set(gca,'Fontsize',fontsize)
colormap(mycolormap);colorbar;
xlabel('Minimum Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
ylabel('(hour$^{-1}$)','interpreter','latex');
title('Growth Rate vs. ${R_i}_\mathrm{min}$ Colored by Shear, $\log(\Lambda)$ (s$^{-1}$)','interpreter','latex');



% Colored by Slope Burger Number
subplot(2,2,4);hold on;
scatter(Ri_sort,grow_sort,markersize,(Bu_sort))
set(gca,'xscale','log')
xlim([0 120])
ylim([0 0.85])
grid on;box on;
set(gca,'Fontsize',fontsize)
colormap(mycolormap);colorbar;
clim([0 3]);
xlabel('Minimum Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
ylabel('(hour$^{-1}$)','interpreter','latex');
title('Growth Rate vs. ${R_i}_\mathrm{min}$ Colored by $B_u=\tilde N\theta/f_0$','interpreter','latex');



print('-dpng','-r300','fig_supp/figS_sort_eig_matlab.png');
