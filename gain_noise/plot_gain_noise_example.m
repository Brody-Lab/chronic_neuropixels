function[] = plot_gain_noise_example(example, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
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
ylabel('RMS noise (uV)')
D.noise_uV(D.noise_uV>P.noise_threshold_uV)=P.noise_threshold_uV;
hdl = plot(D.noise_uV, 'ko', 'MarkerSize', 1);
if example=="unimplanted"
%     set(hdl, 'Color', P.color_order(3,:))
    title_text='Unimplanted';
else
    idx =strcmp(T.recording_id, P.gain_noise_example.(example));
    title_text = ['Explanted after \newline' num2str(T.cumul_days_implanted(idx)) ' days implanted'];
end
title(title_text)
xlabel('deep <- electrodes -> superficial')