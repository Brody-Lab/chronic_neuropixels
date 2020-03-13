% PLOT_GAIN_NOISE_SUMMARY plot the root-mean-square noise of the implanted
% portion of the probe against the number of days the probe has been
% implanted
function[] = plot_gain_noise_summary()
P=get_parameters;
T= readtable(P.gain_noise_log_path);
unique_probes = unique(T.probe_sn);
n_probes = numel(unique_probes);
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
        idx_electrodes = 1:T.electrodes_implanted(ind);
        gain_median{i,1}(j,1) = median(D.gain(idx_electrodes));
        gain_95{i,1}(j,1:2) = quantile(D.gain(idx_electrodes), [0.025, 0.975]);
        noise_uV_median{i,1}(j,1) = median(D.noise_uV(idx_electrodes));
        noise_uV_95{i,1}(j,1:2) = quantile(D.noise_uV(idx_electrodes), [0.025, 0.975]);
        frac_noisy{i,1}(j,1) = sum(D.noise_uV(idx_electrodes)>P.noise_threshold_uV)/numel(idx_electrodes);
    end
end
%% noise
figure('Position', P.figure_position_gn_summary)
set(gca, P.axes_properties_behavior{:}, ...
         'XLim', [-10, 150])
xlabel('Days implanted')
ylabel('Median RMS noise (uV)')
for i =1:n_probes
    plot(days_implanted{i}, noise_uV_median{i}, 'ko--')
end
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.gain_noise_fldr_path filesep 'RMS_noise_summary'], P.figure_image_format{i})
end
%% fraction of electrodes with RMS noise above threshold
figure('Position', P.figure_position_gn_summary)
set(gca, P.axes_properties_behavior{:}, ...
         'XLim', [-10, 150])
xlabel('Days implanted')
ylabel(['Fraction > ' num2str(P.noise_threshold_uV) ' uV'])
for i =1:n_probes
    plot(days_implanted{i}, frac_noisy{i}, 'ko--')
end
set(gca, 'Ylim', ylim.*[0,1])
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.gain_noise_fldr_path filesep 'frac_noisy_electrodes'], P.figure_image_format{i})
end