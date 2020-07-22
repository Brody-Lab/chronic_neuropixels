function[]=plot_longevity_by_brain_area(Cells, brain_area, varargin)
% PLOT_LONGEVITY_BY_BRAIN_AREA Show longevity data from multiple probes for
% a specified brain area.
%
%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILES 
%
%   brain_area
%       A char vector specifying the brain area where recordings are
%       plotted
%
%=OPTIONAL INPUT
%
%   axes
%       An Axes handle where the plot is generated. If it is empty
%       (default), then a new figure is created.
%
%   metric
%       The metric of neuronal signal to be plotted
%
%   varying
%       The feature that is varied for each implant to distinguish it from
%       others. Either the marker type or the color. 
%
%   ylabel_on
%       A scalar logical specifying whether to show the ylabel
%   
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'metric', 'unit', ...
    @(x) validateattributes(x, {'string', 'cell', 'char'}, {}))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'varying', 'marker', @(x) ismember(x, {'color', 'marker'}));
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
brain_area = string(brain_area);
assert(ismember(brain_area, [ "MCtx"
                                "S1"
                                "SC"
                                "amygdala"
                                "dStr"
                                "dlStr"
                                "dmFC"
                                "nIC"
                                "pallidum"
                                "piriform"
                                "postsubiculum"
                                "striatum tail"
                                "vStr"
                                "vmFC"]));
T = get_metrics_from_Cells(Cells, 'condition_on', 'brain_area', ...
                                  'brain_area', {brain_area, {'other'}}, ...
                                  'standardize_area_names', true);
T = T(T.brain_area==brain_area,:);
%% normalize by electrodes
y_label = P.text.(P_in.metric);
if numel(unique(T.condition)) > 1 && ~P_in.normalize_by_electrodes
    P_in.normalize_by_electrodes=true;
    warning('Because of conditioning, normalize_by_electrodes is set to be true')
end
if P_in.normalize_by_electrodes && ~strcmp(P_in.metric, 'Vpp')
    metric=T.(P_in.metric)./T.n_elec;
    y_label = [y_label '/electrode'];
else
    metric=T.(P_in.metric);
end
if P_in.metric == "event_rate"
    y_label = [y_label, ' (Hz)'];
end
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
%%
rats = unique(T.rat)';
for i = 1:numel(rats)
    idx = T.rat == rats(i);
    x = T.days_elapsed(idx);
    y = metric(idx);
    [x, idx] = sort(x);
    y = y(idx);
    hdl(i) = plot(x,y, 'o-');
end
legend(hdl, rats)
ylim([0,1].*ylim)
xlabel('Days since implant')
if P_in.ylabel_on
    ylabel(y_label); 
end