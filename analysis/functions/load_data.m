    
    

clear oceTAUX oceTAUY
%     if(expname(end-4:end)=='_prod')
    if(is_prod_run(ne))
        load([prodir expname '_tavg_5yrs.mat'],'THETA','SALT','UVEL','VVEL','VVELTH','UVELTH','ETAN',...
                'SHIfwFlx','oceTAUX','oceTAUY','SHI_TauX','SHI_TauY','WVEL','WVELTH');
        tt = THETA;
        ss = SALT;
        uu = UVEL;
        vv = VVEL;
        ww = WVEL;
        vt = VVELTH;
        ut = UVELTH;
        wt = WVELTH;
        eta = ETAN;
        %    rho_insitu = RHOAnoma+rhoConst;
        if(useSEAICE)
            load([prodir expname '_tavg_5yrs.mat'],'SIuice','SIvice');
            ui = SIuice;
            vi = SIvice;
        end
    else
        tt = rdmds([exppath,'/results/THETA'],nIter(n));
        ss = rdmds([exppath,'/results/SALT'],nIter(n));
        uu = rdmds([exppath,'/results/UVEL'],nIter(n));
        vv = rdmds([exppath,'/results/VVEL'],nIter(n));
        ww = rdmds([exppath,'/results/WVEL'],nIter(n));
        vt = rdmds([exppath,'/results/VVELTH'],nIter(n));
        ut = rdmds([exppath,'/results/UVELTH'],nIter(n));
        wt = rdmds([exppath,'/results/WVELTH'],nIter(n));
        eta = rdmds([exppath,'/results/ETAN'],nIter(n));
        SHIfwFlx = rdmds([exppath,'/results/SHIfwFlx'],nIter(n));
        SHI_TauY = rdmds([exppath,'/results/SHI_TauY'],nIter(n));
        if(useSEAICE)
            ui = rdmds([exppath,'/results/SIuice'],nIter(n));
            vi = rdmds([exppath,'/results/SIvice'],nIter(n));
        end
    end
    