


period = diff(tt([1 end])) ./ (0:1:Nt/2-1);
freq = 2*pi./period;

energy_detrend=detrend((pe/median(pe)+ke/median(ke))/2);
energy_detrend = smooth(energy_detrend,(43200/dt));

energy_fft = fft(energy_detrend)/length(energy_detrend);
energy_fft = energy_fft/median(energy_fft);

figure(5)
clf
loglog(freq/omega,abs(energy_fft(1:Nt/2)).^2,'LineWidth',2);
hold on;
grid on;grid minor;

%%% Find the dominant low frequency of the energy spectrum
[max_energy, max_energy_idx]=max(abs(energy_fft(1:Nt/2)).^2)
dominant_period = period(max_energy_idx);

Nperiod_dominant = floor(Lt*(1-3/NTtide)/dominant_period)-0.25;
end_idx_tt = floor(dominant_period*Nperiod_dominant/dt);
begin_idx_tt = floor(dominant_period/dt);

    fit_span = begin_idx_tt:end_idx_tt;
    xxplot = tt/3600;
    yyplot = log(pe/median(pe)+ke/median(ke))/2;
    yyplot = smooth(yyplot,(43200/dt));
    % yyplot = log(pe+ke)/2;
    [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    grow(i) = pp(1)
    if(isnan(grow(i)))
        warning('NaN in growth rate!')
    end
    [y_fit,delta_fit] = polyval(pp,xxplot,S);
    figure(20)
    clf;
    plot(xxplot,yyplot)
    hold on;grid on;grid minor;
    plot(xxplot(fit_span), y_fit(fit_span));
    scatter(xxplot(begin_idx_tt),y_fit(begin_idx_tt));
    scatter(xxplot(end_idx_tt), y_fit(end_idx_tt));
    hold off;