    %%% Start the loop
function [grow,bb,uu,vv,ww,pp,re_bb,re_uu,re_vv,re_ww,ct,st,mz_t,a1_t,ke_nond,grav,pe_nond,fit_span,xxplot,yyplot]...
        =loop(grow,j,NTtide,dt,Nt,dbdt,dudt,dvdt,dwdt,omega,m0,k0,shear,ss,cs,N,kappa,nu,tt,bb,uu)

    for o=1:Nt-1

        %%% Fourth-order Runge-Kutta method %%%
        RKstep = 1;
        t0 = tt(o);
        u0 = uu(o);
        v0 = vv(o);
        b0 = bb(o);
        [ww,pp,dwdt,dudt,dvdt,dbdt]=tendency(o,RKstep,ww,pp,dwdt,dudt,dvdt,dbdt,t0,m0,k0,l0,u0,v0,b0,ss,cs,shear,omega,N,f0,kappa,nu);
        k_1u = dudt(o);
        k_1v = dvdt(o);
        k_1b = dbdt(o);

        % Euler forward predictor advancing dt/2:
        RKstep = 2;
        b_2 = bb(o)+0.5*dt*k_1b;
        u_2 = uu(o)+0.5*dt*k_1u;
        v_2 = vv(o)+0.5*dt*k_1v;
        t0 = tt(o)+dt/2;
        b0 = b_2;
        u0 = u_2;
        v0 = v_2;
        [ww,pp,dwdt,dudt,dvdt,dbdt]=tendency(o,RKstep,ww,pp,dwdt,dudt,dvdt,dbdt,t0,m0,k0,l0,u0,v0,b0,ss,cs,shear,omega,N,f0,kappa,nu);        k_2w = dwdt(o);
        k_2u = dudt(o);
        k_2v = dvdt(o);
        k_2b = dbdt(o);

        % Euler backward corrector advancing dt/2:
        RKstep = 3;
        b_3 = bb(o)+0.5*dt*k_2b;
        u_3 = uu(o)+0.5*dt*k_2u;
        v_3 = vv(o)+0.5*dt*k_2v;
        t0 = tt(o)+dt/2;
        b0 = b_3;
        u0 = u_3;
        v0 = v_3;
        [ww,pp,dwdt,dudt,dvdt,dbdt]=tendency(o,RKstep,ww,pp,dwdt,dudt,dvdt,dbdt,t0,m0,k0,l0,u0,v0,b0,ss,cs,shear,omega,N,f0,kappa,nu);        k_3w = dwdt(o);
        k_3b = dbdt(o);
        k_3u = dudt(o);
        k_3v = dvdt(o);
        
        % Mid-point predictor advancing dt:
        RKstep = 4;
        b_4 = bb(o)+dt*k_3b;
        u_4 = uu(o)+dt*k_3u;
        v_4 = vv(o)+dt*k_3v;
        t0 = tt(o)+dt;
        b0 = b_4;
        u0 = u_4;
        v0 = v_4;
        [ww,pp,dwdt,dudt,dvdt,dbdt]=tendency(o,RKstep,ww,pp,dwdt,dudt,dvdt,dbdt,t0,m0,k0,l0,u0,v0,b0,ss,cs,shear,omega,N,f0,kappa,nu);        k_4b = dbdt(o);
        k_4b = dbdt(o);
        k_4u = dudt(o);
        k_4v = dvdt(o);
    
        % Simpson rule corrector advancing dt:
        bb(o+1) = bb(o) + (1/6)*(k_1b+2*k_2b+2*k_3b+k_4b)*dt;
        uu(o+1) = uu(o) + (1/6)*(k_1u+2*k_2u+2*k_3u+k_4u)*dt;
        vv(o+1) = vv(o) + (1/6)*(k_1v+2*k_2v+2*k_3v+k_4v)*dt;

    end

   

    re_bb = real(bb);
    re_uu = real(uu);
    re_ww = real(ww);
    re_vv = real(vv);
    re_pp = real(pp);

    pe = re_bb.^2;
    ke = 0.5*(re_uu.^2+re_vv.^2+re_ww.^2);
    % %%% To match Radko (2019) Eq.(19)

    if(nu*kappa~=0)
        Pr = nu/kappa;
    else
        Pr = 7;
    end

    % ke_nond = (k0^2+mz_t.^2).*abs(psi/kappa_const).^2/4;    %%% To match Radko (2019) Eq.(19) 
    % grav = 10;
    % pe_nond = Pr*(abs(bb)/grav).^2/4; %%% Non-dimensionalized KE and PE

    fit_span = Nt/NTtide*1+1:Nt;
    xxplot = tt/3600;
    % yyplot = log(ke_nond+pe_nond)/2;
    yyplot = log(ke+Pr*pe)/2;
    [gr,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    grow(j)=gr(1);
    % gr(1)
    % if(isnan(gr(1)))
    %     warning('NaN in growth rate!')
    % end
    % 
    % [y_fit,delta_fit] = polyval(gr,xxplot,S);
    % fig = figure(20);
    % clf;set(gcf,'Color','w')
    % plot(xxplot/24,yyplot,'LineWidth',2)
    % hold on;grid on;grid minor;
    % plot(xxplot(fit_span)/24, y_fit(fit_span));
    % hold off;
    % xlabel('Time (days)')
    % set(gca,'Fontsize',20)

    % saveas(fig,[[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '.jpeg']);


end