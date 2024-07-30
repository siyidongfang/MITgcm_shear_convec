
clear;close all

expdir = 'backup_exps/parallel_flat_nu2e-6/';
shear_all = (0:0.1:2.0)/1e3;

h_shear = 1000;
m0_limit = 2*pi/h_shear;

m0max = 2*pi/1;
m0min = m0_limit;
k0max = 2*pi/3;
k0min = 2*pi/30000;
m0_all = [m0min:0.01/10:0.6 0.61:0.01/2:1];
kx_all = [k0min:0.001/10:0.06 0.061:0.001/2:0.1];
% k0min = 2*pi/100000;
% m0_all = [m0min 0.026:0.01/10:0.6 0.61:0.01/2:1];
% kx_all = [k0min*[1:1/2:15] 0.001:0.001/10:0.06 0.061:0.001/2:0.1];
Nk = length(kx_all);
Nm = length(m0_all);
lam_z_all = 2*pi./m0_all;
lam_x_all = 2*pi./kx_all;

mmidx = 1:Nm;
kkidx = 20:Nk;
% kkidx = 16:Nk;

omega = 2*pi/43200;

grow_skm = NaN*zeros(length(shear_all),length(kx_all),length(m0_all));

h_shear = 250;
m0_limit = 2*pi/h_shear;

for s = 1:length(shear_all)
% for s = 8
    shear = shear_all(s)
    for i=kkidx
        kx = kx_all(i);
        fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_kx' num2str(i) '.mat'];
        if(isfile(fname))
            load(fname,'grow');
            m0_limit_t = kx *shear/omega + m0_limit;
            % m0_limit_t = 0;
    
            for j=1:Nm
                m0 = m0_all(j);
                if(m0>=m0_limit_t)
                    grow_skm(s,i,j) = grow(j);
                end
            end
        end

    end


    figure(1)
    set(gcf,'Color','w')
    pcolor(kx_all,m0_all,squeeze(grow_skm(s,:,:))');shading flat;colorbar;colormap(redblue);
    clim([-0.3 0.3])
    set(gca,'FontSize',20);ylabel('k (1/m)');xlabel('m (1/m)')
    title('Growth rate (1/hour)')


    max_grow(s) = max(grow_skm(s,:,:),[],'all');

end

figure(2)
set(gcf,'Color','w')
hold on;
scatter(shear_all,max_grow);

shear_km = shear_all;
growth_km = max_grow;

%%% Find the wavenumber of the most unstale mode

mm_crop = m0_all(mmidx);
kk_crop = kx_all(kkidx);

for s = 1:length(shear_km)
   aaa = squeeze(grow_skm(s,kkidx,mmidx));


    [growth_crop(s) I] = max(aaa,[],'all');
    [kk_idx(s),mm_idx(s)] = ind2sub(size(aaa),I);
    max_m0(s) = mm_crop(mm_idx(s));
    max_kx(s) = kk_crop(kk_idx(s));
end

max_lambda_z = 2*pi./max_m0;
max_lambda_x = 2*pi./max_kx;

min_mt = (max_kx.*shear_all/omega + m0_limit);

% save('growth_experiments_flat_nu2e-4.mat','shear_km','growth_km','kx_all','m0_all','grow_smk',...
%     'mmidx','kk_idx','max_lambda_z','max_lambda_x','max_m0','max_kx')

max_m0(isnan(max_grow))=NaN;
max_kx(isnan(max_grow))=NaN;

% save('parallel_flat_nu2e-6.mat','max_m0','max_kx','max_grow');
