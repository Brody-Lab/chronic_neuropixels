% PLOT_GAIN_NOISE_BROKEN_FRAC The fraction of electrodes that are above a
% noise threshold
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = plot_gain_noise_broken_frac(varargin)
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
%% fraction of electrodes with RMS noise above threshold
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
    hdl=plot(days_implanted{i}, frac_noisy{i}, 'ko--', 'linewidth', 1);
    if days_implanted{i} == 0
        set(hdl, 'marker', '*');
    end
end
set(gca, 'Ylim', ylim.*[0,1])
xlabel('Cumulative days implanted')
ylabel(['Fraction > ' num2str(P.noise_threshold_uV) ' uV'])
%% stats
idx_new = cell2mat(days_implanted) == 0;
noise_uV_all=vertcat(noise_uV{:});
noise_uV_new = cell2mat(noise_uV_all(idx_new));
noise_uV_exp = cell2mat(noise_uV_all(~idx_new));
x.new(1,1) = sum(noise_uV_new<P.noise_threshold_uV);
x.new(2,1) = sum(noise_uV_new>=P.noise_threshold_uV);
x.exp(1,1) = sum(noise_uV_exp<P.noise_threshold_uV);
x.exp(2,1) = sum(noise_uV_exp>=P.noise_threshold_uV);
x=struct2table(x);
[~,pval]=fishertest(x);
[~,ci_new] = binofit(x.new(2), sum(x.new));
[~,ci_exp] = binofit(x.exp(2), sum(x.exp));
fprintf('\nComparing the fraction of electrodes with a RMS noise > 20 uV between new and explanted probes:')
fprintf('\n   new: %0.2f,95%%CI: [%0.2f, %0.2f]', x.new(2)/sum(x.new)*100, ci_new(1)*100, ci_new(2)*100)
fprintf('\n   exp: %0.2f 95%%CI: [%0.2f, %0.2f]', x.exp(2)/sum(x.exp)*100, ci_exp(1)*100, ci_exp(2)*100)
fprintf('\n   (using the Fisher'' exact test):')
fprintf('\n   p = %f', pval)
