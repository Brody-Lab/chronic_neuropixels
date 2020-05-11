




analyze_gain_noise_data;
all_noise_new=cellfun(@(x)x(end),noise_uV(idx_new));
all_noise_new = cat(1,all_noise_new{:});
all_noise_exp=cellfun(@(x)x(end),noise_uV(~idx_new));
all_noise_exp = cat(1,all_noise_exp{:});

figure;
lims=[5 20];
factor=1.02;
all_noise_new = max(all_noise_new,lims(1));
all_noise_new = min(all_noise_new,lims(2));
all_noise_exp = max(all_noise_exp,lims(1));
all_noise_exp = min(all_noise_exp,lims(2));

histogram_log((all_noise_new),'NumBins',500,'Normalization','cdf','DisplayStyle','stairs','BinLimits',[1/factor factor].*lims);
hold on;histogram_log((all_noise_exp),'NumBins',500,'Normalization','cdf','DisplayStyle','stairs','BinLimits',[1/factor factor].*lims);
set(gca,'xlim',[1/factor factor].*lims,'ylim',[0 1],'ygrid','on','xgrid','on');


