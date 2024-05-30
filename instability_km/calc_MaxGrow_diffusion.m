
% clear;close all

expdir = 'experiments_flat_nu2e-4/';
shear_all = (0:0.1:2.0)/1e3;

m0max = 2*pi/1;
m0min = 2*pi/4000;
k0max = 2*pi/3;
k0min = 2*pi/100000;
m0_all = [0 m0min m0min*2 m0min*3 m0min*4 0.01:0.01:1];
kx_all = [0 k0min k0min*2 k0min*3 k0min*4 0.001:0.001:0.1];

lam_z_all = 2*pi./m0_all;
lam_x_all = 2*pi./kx_all;

mmidx = 8:105;
kkidx = 7:105;

omega = 2*pi/43200;



grow_smk = NaN*zeros(length(shear_all),length(m0_all),length(kx_all));
grow_smk_limit = NaN*zeros(length(shear_all),length(m0_all),length(kx_all));

h_shear = 250;
m0_limit = 2*pi/h_shear;

for s = 1:length(shear_all)
    shear = shear_all(s)
    for m=1:length(m0_all)
        m0 = m0_all(m);
        fname = [[expdir 'shear_' num2str(shear*1e3,3)] '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.mat'];
        if(isfile(fname))
            load(fname,'grow');
            grow_smk(s,m,:) = grow;
            grow_smk_limit(s,m,:) = grow;
        end

        if(m0<m0_limit)
            grow_smk_limit(s,m,:) = NaN;
        end

        for k=1:length(kx_all)
            kx=kx_all(k);
            m0_limit_t = kx *shear/omega + m0_limit;
            if(m0<m0_limit_t)
                grow_smk_limit(s,m,k) = NaN;
            end
        end

    end


    figure(1)
    set(gcf,'Color','w')
    pcolor(kx_all(kkidx),m0_all(mmidx),squeeze(grow_smk_limit(s,mmidx,kkidx)));shading flat;colorbar;colormap(redblue);
    clim([-0.4 0.4])
    set(gca,'FontSize',20);xlabel('k (1/m)');ylabel('m (1/m)')
    title('Growth rate (1/hour)')
    % ylim([0 1]);xlim([0 0.3]);


    max_grow(s) = max(grow_smk_limit(s,mmidx,kkidx),[],'all');

end

figure(2)
set(gcf,'Color','w')
hold on;
plot(shear_all,max_grow);

shear_km = shear_all;
growth_km = max_grow;


for s = 1:length(shear_km)
   aaa = squeeze(grow_smk_limit(s,mmidx,kkidx));
   mm_crop = m0_all(mmidx);
   kk_crop = kx_all(kkidx);

   [growth_crop(s) I] = max(aaa,[],'all');
   [mm_idx(s),kk_idx(s)] = ind2sub(size(aaa),I);
   max_m0(s) = mm_crop(mm_idx(s));
   max_kx(s) = kk_crop(kk_idx(s));
end

max_lambda_z = 2*pi./max_m0;
max_lambda_x = 2*pi./max_kx;

figure(10)
hold on;
plot(shear_all,growth_crop);

% save('growth_experiments_flat_nu2e-4.mat','shear_km','growth_km','kx_all','m0_all','grow_smk',...
%     'mmidx','kk_idx','max_lambda_z','max_lambda_x','max_m0','max_kx')

