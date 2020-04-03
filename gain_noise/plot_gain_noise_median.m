% PLOT_GAIN_NOISE_MEDIAN the median noise level
function[] = plot_gain_noise_median(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
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
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:}, ...
         'XLim', [-10, 150])
for i =1:n_probes
    hdl = plot(days_implanted{i}, noise_uV_median{i}, 'ko--');
    if days_implanted{i} == 0
        set(hdl, 'marker', '*');
    end
end
xlabel('Cumulative days implanted')
ylabel('Median RMS noise (uV)')