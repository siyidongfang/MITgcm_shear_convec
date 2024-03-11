function make_TS_plot_gsw (salt,pot_temp,depths,p_ref,pot_dens_contours,lat_ref,lon_ref)
%%%
%%% USAGE: make_TS_plot_gsw (salt, pot_temp, depths, p_ref, pot_dens_contours)
%%%
%%% Makes a temperature/salinity plot using the 2D salinity/potential 
%%% temperature data, using the GSW Toolbox.
%%%
%%% Arguments:
%%% salt - salinity matrix
%%% pot_temp - potential temperature matrix
%%% depths - vector of depths
%%% p_ref - scalar reference pressure for contouring density lines
%%% pot_dens_contours - vector of potential density lines to plot.
%%%                     If a single integer, e.g. 30, is specified then
%%%                     in this example 30 evenly-spaced contours will be plotted.
%%%    
%%% salt and pot_temp must be the same size (M x N matrices). depths must 
%%% be a vector of length N. p_ref must be a positive pressure in decibars.
%%%
%%% lon_ref =  longitude in decimal degrees                      [ 0 ... +360 ]
%%%                                                   or  [ -180 ... +180 ]
%%% lat_ref  =  latitude in decimal degrees north                [ -90 ... +90 ] 
%%% lat_ref & lon_ref must have the dimension of 1x1.


%%% Error checking
if (size(salt) ~= size(pot_temp))
  error('make_TS_plot: Sizes of salinity and potential temperature input data must be the same');
end
if (length(depths) ~= size(salt,2))
  error('make_TS_plot: Length of depth data must match number of columns of salinity and temperature data');
end
if (p_ref < 0)
  error('make_TS_plot: Reference pressure must be >= 0');
end

%%% Easier variable names
pt = pot_temp;
ss = salt;
Ny = size(ss,1);
Nz = size(ss,2);

%%% Get ranges of T/S
pt_max = max(max(pt)) + 1;
pt_min = min(min(pt)) - 1;
ss_max = max(max(ss)) + 0.1;
ss_min = min(min(ss)) - 0.1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Calculate potential density
% pd = densmdjwf(SS_grid,PT_grid,p_ref) - 1000;
%%% Use the GSW toolbox instead of densmdjwf
SA_grid = gsw_SA_from_SP(SS_grid,p_ref,lon_ref,lat_ref);  
CT_grid = gsw_CT_from_pt(SA_grid,PT_grid); 
pd = gsw_rho(SA_grid,CT_grid,p_ref) - 1000;

%%% Plotting options
scrsz = get(0,'ScreenSize');
fontsize = 22;
plotloc = [0.1 0.15 0.8 0.8];
framepos = [scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2];

%%% Scatter the T/S data
handle = figure;
set(handle,'Position',framepos);
clf;
for j=1:Ny
  scatter(ss(j,:),pt(j,:),4,-depths);
  if (j==1)
    hold on;
  end
end
[C,h] = contour(SS_grid,PT_grid,pd,pot_dens_contours,'EdgeColor','k');
clabel(C,h);
% hold off;
xlabel('Salinity (g/kg)','interpreter','latex');
ylabel('Potential temperature ($^\circ$C)','interpreter','latex');
set(gca,'Position',plotloc);
handle = colorbar;
colormap(flipdim(jet(100),2));
set(handle,'FontSize',fontsize);
caxis([min(-depths) max(-depths)]);
%axis([ss_min ss_max pt_min pt_max]);
% axis([33.4 34.7 -2 -1]);

annotation('textbox',[0.82 0.05 0.3 0.05],'String','Depth (m)','interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
set(gca,'FontSize',fontsize);
