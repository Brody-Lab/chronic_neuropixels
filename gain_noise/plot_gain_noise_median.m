% PLOT_GAIN_NOISE_MEDIAN plot the median the median noise level
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
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
        noise_uV{i,1}{j,1} = D.noise_uV(idx_electrodes);
        noise_uV_median{i,1}(j,1) = median(D.noise_uV(idx_electrodes));
        noise_uV_95{i,1}(j,1:2) = quantile(D.noise_uV(idx_electrodes), [0.025, 0.975]);
        frac_noisy{i,1}(j,1) = sum(D.noise_uV(idx_electrodes)>P.noise_threshold_uV)/numel(idx_electrodes);
        days_implanted{i}(days_implanted{i}>200) = 200;
    end
end
%% plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:}, ...
         'XLim', [-10, 200], ...
         'Xtick', [0, 100, 200], ...
         'XTicklabel', {'0', '100', '\geq200'})
for i =1:n_probes
    hdl = plot(days_implanted{i}, noise_uV_median{i}, 'ko--', 'linewidth', 1);
    if days_implanted{i} == 0
        set(hdl, 'marker', '*');
        hdl_new = hdl;
    else
        hdl_explanted = hdl;
    end
%     if numel(noise_uV_median{i}) > 1
%         for j = 1:numel(noise_uV_median{i})
%             plot(days_implanted{i}(j), noise_uV_median{i}(j), 'o', 'Color', P.color_order(j,:), 'linewidth', 1);
%         end
%     end
end
xlabel('Cumulative days implanted')
ylabel('Median noise (uV_{RMS})')
% legend([hdl_new, hdl_explanted], {'New', 'Explanted'}, 'Location', 'Best')
%% stats: Compare between fresh and explanted probes
fprintf('\nComparing the RMS noise on the electrodes between the new and explanated probes:')
days_implanted_latest = cellfun(@(x) x(end), days_implanted);
idx_new = days_implanted_latest == 0;
noise_uV_all=cellfun(@(x) x(end), noise_uV);
noise_uV_new = cell2mat(noise_uV_all(idx_new));
noise_uV_exp = cell2mat(noise_uV_all(~idx_new));
fprintf('\n   Median=%0.3f (%i new probes, n = %i electrodes) and %0.3f (%i explanated, n = %i electrodes)', ...
            median(noise_uV_new), ...
            sum(idx_new), ...
            numel(noise_uV_new), ...
            median(noise_uV_exp), ...
            sum(~idx_new), ...
            numel(noise_uV_exp))
pval=ranksum(noise_uV_new, noise_uV_exp);
fprintf('\n   p = %e\n', pval)
%% stats; For each probe that was implanted twice, do some statistics
idx = unique_probes == 17131312432;
fprintf('\nComparing the RMS noise one one probe across multiple explantations:')
for i = 1:numel(noise_uV{idx})
    fprintf('\n   After %i explantation, the median noise = %0.3f uV.', i, median(noise_uV{idx}{i}))
end
% do an anova
% make the dependent variable
d=[];
for i = 1:numel(noise_uV{idx})
    d=[d;repmat(days_implanted{idx}(i), numel(noise_uV{idx}{i}), 1)];
end
pval = anovan(cell2mat(noise_uV{idx}), d, 'disp', 'off');
fprintf('\n   ANOVA on RMS noise with cumul. days implanted as a main factor:')
fprintf('\n      p = %0.1e', pval)
fprintf('\n')
%% wost vs best
median_new = cellfun(@median, noise_uV_all(idx_new));
median_exp = cellfun(@median, noise_uV_all(~idx_new));
fprintf('\nDifference in median uV between the best unimplanted and worst explanted probe:')
fprintf('\n    %0.2f     uV', max(median_exp)-min(median_new))