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
%   T_implants
%       If this were provided, then the marker/color specs are provided in
%       the order of the implants table
%
%   metric
%       The metric of neuronal signal to be plotted
%
%   time_bin_edges
%       Bin edges for discretize data
%
%   varying
%       The feature that is varied for each implant to distinguish it from
%       others. Either the marker type or the color. 
%
%   legend_on
%       A scalar logical specifying whether to show the leged
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
addParameter(parseobj, 'time_bin_edges', 2.^(-0.5:9.5), ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}))
addParameter(parseobj, 'varying', 'marker', @(x) ismember(x, {'color', 'marker'}));
addParameter(parseobj, 'legend_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
brain_area = string(brain_area);
assert(ismember(brain_area, P.brain_areas));
T_sess = get_metrics_from_Cells(Cells, 'condition_on', 'brain_area', ...
                                      'brain_area', {brain_area, {'other'}}, ...
                                      'standardize_area_names', true);
T_sess = T_sess(T_sess.brain_area==brain_area,:);
%% normalize by electrodes
y_label = P.text.(P_in.metric);
if numel(unique(T_sess.condition)) > 1 && ~P_in.normalize_by_electrodes
    P_in.normalize_by_electrodes=true;
    warning('Because of conditioning, normalize_by_electrodes is set to be true')
end
if P_in.normalize_by_electrodes && ~strcmp(P_in.metric, 'Vpp')
    metric=T_sess.(P_in.metric)./T_sess.n_elec;
    y_label = [y_label '/electrode'];
else
    metric=T_sess.(P_in.metric);
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
set(gca, 'XLim', [min(P_in.time_bin_edges), max(P_in.time_bin_edges)], ...
         'XTick', 2.^(0:9), ...
         'YLim', [0,1.5], ...
         'YTick', 0:0.5:1.5)
%%
T_implants = readtable(P.implants_by_area_path);
T_implants = T_implants(T_implants.brain_area == brain_area,:);
for i = 1:size(T_implants,1)
    idx = strcmp(T_sess.rat, T_implants.rat(i)) & ...
          T_sess.probe_serial == T_implants.probe_serial(i);
    x = groupsummary(T_sess.days_elapsed(idx), T_sess.days_elapsed(idx), ...
                     P_in.time_bin_edges, 'mean');
    y = groupsummary(metric(idx), T_sess.days_elapsed(idx), ...
                     P_in.time_bin_edges, 'mean');    
    [marker, colr] = lookup_linespec_of_implant(T_implants.rat(i), T_implants.probe_serial(i));
    hdl(i) = plot(x,y, marker, 'Color', colr, 'linewidth', 1.5, 'markersize', 10);
end
if P_in.legend_on
    legend(hdl, num2str(T_implants.implant_number), 'location', 'best')
end
ylim([0,1].*ylim)
xlabel('Days since implant')
if P_in.ylabel_on
    ylabel(y_label); 
end