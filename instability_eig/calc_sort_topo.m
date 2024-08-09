clear all;close all;

load('grow_K1.mat')
% load('grow_M2.mat')
% mycolormap = WhiteBlueGreenYellowRed(9);
mycolormap = WhiteBlueGreenYellowRed(8);
fontsize = 20;markersize = 36;

figure(1);clf;
set(gcf,'Color','w')
set(gcf,'Color','w','Position',[64 68 1200 950]);

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

    subplot(3,1,1);hold on;
    scatter(Ri_sort,grow_sort,markersize,topo.*ones(1,length(grow_sort)),'filled')
    % plot(Ri_sort,grow_sort)
    set(gca,'xscale','log')
    xlim([0 120])
    ylim([0 0.85])
    grid on;box on;
    set(gca,'Fontsize',fontsize)
    colormap(mycolormap);colorbar;
    clim([0 20]);
    xlabel('Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
    ylabel('(hour$^{-1}$)','interpreter','latex');

    title('Growth Rate vs. ${R_i}_\mathrm{min}$ Colored by Topographic Slope ($^\circ$)','interpreter','latex');

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
    
    % subplot(2,2,2);hold on;
    % scatter(1./Ri_sort,grow_sort,markersize,N.*ones(1,length(grow_sort)),'filled')
    % xlim([0 10])
    % ylim([0 0.85])
    % grid on;box on;
    % set(gca,'Fontsize',fontsize)
    % colormap(mycolormap);colorbar;
    % % clim([0 20]);
    % xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
    % ylabel('(hour$^{-1}$)','interpreter','latex');

    subplot(3,1,2);hold on;
    scatter(Ri_sort,grow_sort,markersize,N.*ones(1,length(grow_sort)),'filled')
    set(gca,'xscale','log')
    xlim([0 120])
    ylim([0 0.85])
    grid on;box on;
    set(gca,'Fontsize',fontsize)
    colormap(mycolormap);colorbar;
    % clim([0 20]);
    xlabel('Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
    ylabel('(hour$^{-1}$)','interpreter','latex');
    title('Growth Rate vs. ${R_i}_\mathrm{min}$ Colored by Stratification (s$^{-1}$)','interpreter','latex');

end


% Colored by Shear



