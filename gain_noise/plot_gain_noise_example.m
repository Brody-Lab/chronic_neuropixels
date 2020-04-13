% PLOT_GAIN_NOISE_EXAMPLE plot the input-referred noise of an example probe
%
%=INPUT
%
%   example
%       A char array specifying whether an "implanted" or "unimplanted"
%       probe is shown
%
%=OPTIONAL INPUT
%
%   axes
%       The axes object where the plot is made
%
%   ylabel
%       A logical scalar indicating whether to show the y-axis label

function[] = plot_gain_noise_example(example, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'ylabel', true, @(x) islogical(x) && isscalar(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
assert(any(strcmp(example, {'implanted', 'unimplanted'})))
P=get_parameters;
data_file_path = [P.gain_noise_fldr_path filesep P.gain_noise_example.(example) '.csv'];
D = readtable(data_file_path);
T = readtable(P.gain_noise_log_path);
yticklabel = cellfun(@num2str, num2cell(0:5:P.noise_threshold_uV), 'uni', 0);
yticklabel{end} = ['\geq' num2str(P.noise_threshold_uV)];
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, 'Xlim', [0 961], ...
         'xtick', [], ...
         'YLim', [0, P.noise_threshold_uV], ...
         'Ytick', 0:5:P.noise_threshold_uV, ...
         'YtickLabel', yticklabel)
if P_in.ylabel
    ylabel('Input-referred noise (uV_{RMS})')
end
D.noise_uV(D.noise_uV>P.noise_threshold_uV)=P.noise_threshold_uV;
hdl = plot(D.noise_uV, 'ko', 'MarkerSize', 1);
if example=="unimplanted"
    title_text='Unimplanted';
    fprintf('\nExample probe that has never been implanted')
    fprintf('\n   Median = %0.1f (uV) and CI = [%0.1f, %0.1f] (uV)\n', ...
            median(D.noise_uV), quantile(D.noise_uV, 0.025), quantile(D.noise_uV, 0.975));
else
    idx =strcmp(T.recording_id, P.gain_noise_example.(example));
    title_text = ['Explanted after \newline' num2str(T.cumul_days_implanted(idx)) ' days implanted'];
    fprintf('\nExample probe explanted after %i days implanted', T.cumul_days_implanted(idx))
    fprintf('\n   Median = %0.1f (uV) and CI = [%0.1f, %0.1f] (uV)\n', ...
            median(D.noise_uV), quantile(D.noise_uV, 0.025), quantile(D.noise_uV, 0.975));
end
title(title_text)
xlabel('deep <- electrodes -> superficial')
