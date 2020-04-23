if ~exist('max_z','var')
    max_z=Inf;
end
P=get_parameters;
T= readtable(P.gain_noise_log_path);
unique_probes = unique(T.probe_sn);
n_probes = numel(unique_probes);
min_channels=50;
%max_z=4;
days_implanted = cell(n_probes,1);
T.cumul_days_implanted(isnan(T.cumul_days_implanted)) = 0;
for i = 1:n_probes
    idx = T.probe_sn == unique_probes(i);
    days_implanted{i} = T.cumul_days_implanted(idx);
    inds = find(idx);
    for j = 1:numel(inds)
        ind = inds(j);
        data_file_path = [P.gain_noise_fldr_path filesep T.recording_id{ind} '.csv'];
        D = readtable(data_file_path);
        if i==5
            T.electrodes_implanted(ind) = T.electrodes_implanted(ind) - 10;
        end
        idx_electrodes = 1:T.electrodes_implanted(ind);
        idx_channels = ((1:384) + 384)<=T.electrodes_implanted(ind);
        if sum(idx_channels)<min_channels
            bank_corr_rank{i,1}(j,1) = NaN;
            bank_corr_linear{i,1}(j,1) = NaN;       
            bank_corr_rank_z{i,1}(j,1) = NaN;
            bank_corr_linear_z{i,1}(j,1) = NaN;              
            vare_mean{i,1}(j,1) = NaN;    
            bank_0_noise{i,1} = NaN;
            bank_1_noise{i,1} = NaN;
            bank_0_noise_z{i,1} = NaN;
            bank_1_noise_z{i,1} = NaN;  
            bank_rsquare_z{i,1} = NaN;            
            continue
        end
        bank_idx =1:384;
        bank_0_noise{i,1} = table2array(D(bank_idx(idx_channels),2));
        bank_1_noise{i,1} = table2array(D(bank_idx(idx_channels)+384,2));
        bank_0_noise_z{i,1} = zscore(bank_0_noise{i,1});
        bank_1_noise_z{i,1} = zscore(bank_1_noise{i,1});
        outliers = abs(bank_0_noise_z{i,1})>max_z | abs(bank_1_noise_z{i,1})>max_z;
        bank_0_noise{i,1}(outliers) = [];
        bank_1_noise{i,1}(outliers) = [];         
        bank_0_noise_z{i,1}(outliers) = [];
        bank_1_noise_z{i,1}(outliers) = []; 
        %bank_0_noise_z{i,1} = zscore(bank_0_noise_z{i,1});
        %bank_1_noise_z{i,1} = zscore(bank_1_noise_z{i,1});
        bank_corr_rank{i,1}(j,1) = corr(bank_0_noise{i,1},bank_1_noise{i,1},'Type','Spearman');
        bank_corr_linear{i,1}(j,1) = corr(bank_0_noise{i,1},bank_1_noise{i,1});
        mean_bank_noise = (bank_0_noise{i,1}+bank_1_noise{i,1})/2;
        vare_mean{i,1}(j,1) = rsquare(zscore(bank_0_noise{i,1}),zscore(bank_1_noise{i,1}));
        days_implanted{i}(days_implanted{i}>200) = 200;
        bank_corr_rank_z{i,1}(j,1) = corr(zscore(bank_0_noise{i,1}),zscore(bank_1_noise{i,1}),'Type','Spearman');
        bank_corr_linear_z{i,1}(j,1) = corr(zscore(bank_0_noise{i,1}),zscore(bank_1_noise{i,1}));      
        bank_rsquare_z{i,1}(j,1) = rsquare(zscore(bank_0_noise{i,1}),zscore(bank_1_noise{i,1}));        
        
        probe_sn{i,1} = T.probe_sn(ind);
    end
end

days_implanted_latest = cellfun(@(x) x(end), days_implanted);
idx_new = days_implanted_latest == 0;