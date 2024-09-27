clear all;
load('../instability_km/budget/topo4_output_R3.mat')
% load('../instability_km/budget/flat_output_R3.mat')

Nall = length(re_uuu);
Nend = Nall/140*2;
u = re_uuu(1:Nend);
w = re_www(1:Nend);
t = tt(1:Nend);
u_int = cumsum(u.*t);
w_int = cumsum(w.*t);

t_plot = tt(1:Nall/140);
t_plot = [t_plot t_plot]

figure(1);clf;
plot(t,w_int)
hold on;
plot(t,u_int)
grid on;grid minor;


figure(2);clf;
scatter(u_int,w_int,20,t_plot/43200)
grid on;grid minor;
colorbar
