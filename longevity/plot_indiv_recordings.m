% PLOT_INDIVI_RECORDINGS Show the individul recordins from one of three
% brain regions: 1) medial frontal cortex, 2) motor cortex and dorsal
% striatum, 3) ventral striatum
%
%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILES 
%
%   brain_region
%       A cell, string, or char array specifying the brain region from
%       which the recordings are shown. The options are 'mFC', 'MCtx_ADS',
%       or 'AVS'.
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
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
function [] = plot_indiv_recordings(Cells, brain_area, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'metric', 'single_unit', ...
    @(x) validateattributes(x, {'string', 'cell', 'char'}, {}))
addParameter(parseobj, 'varying', 'marker', @(x) ismember(x, {'color', 'marker'}));
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'legend_on', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
assert(ismember(brain_area, {'vmFC', 'MCtx', 'vmStr'}))
switch brain_area
    case 'vmFC'
        rat_name = ["T212"; "T224"; "T249"];
        unique_bank = 0;
        title_text = 'Ventromedial frontal ctx';
    case 'ADS'
        rat_name = ["T181"; "T182"; "T219"];
        unique_bank = 1;
        title_text = 'Dorsal striatum';
    case 'MCtx_ADS'
        rat_name = ["T181"; "T182"; "T219"];
        unique_bank = 1;
        title_text = 'Motor ctx and dorsal striatum';
    case 'MCtx'
        rat_name = ["T181"; "T182"; "T219"];
        unique_bank = 1;
        title_text = 'Motor cortex';
    case 'vmStr'
        rat_name = ["T181"; "T182"; "T219"];
        unique_bank = 0;
        title_text = 'Ventromedial striatum';
end
P = get_parameters;
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
Cells = standardize_brain_area_names(Cells);
for i = 1:numel(rat_name)
    if ~isempty(brain_area)
        T = get_metrics_from_Cells(Cells, 'condition_on', {'EI', 'brain_area'}, ...
                                          'brain_area', {{brain_area}, {'other'}}, ...
                                          'EI_bin_edges', unique_bank*384 + [1, 384]);
        T = T(contains(T.brain_area, brain_area),:);
    else
        T = get_metrics_from_Cells(Cells, 'condition_on', 'EI', ...
                                          'EI_bin_edges', unique_bank*384 + [1, 384]);
    end
    T = T(T.rat==rat_name(i),:);
    T = T(T.n_elec > 1,:);
    % if there are multiple probes per animal, use only the one with the
    % largest number of recordings
    Ct = groupcounts(T, 'probe_serial'); 
    T = T(T.probe_serial == Ct.probe_serial(end), :);
    switch P_in.varying
        case 'color'
            hdl(i)= plot(T.days_elapsed, T.(P_in.metric)./T.n_elec, ...
                 'o-', 'Color', P.color_order(i,:));
        case 'marker'
            hdl(i) = plot(T.days_elapsed, T.(P_in.metric)./T.n_elec, ...
                 ['k', P.marker_order(i), P.line_order{i}], 'linewidth', 1);
        otherwise
            error('This feature has not been implemented for distinguishing among recordings')
    end 
end
set(gca, 'xlim', [2^-0.5, 2^9.1], ...
         'xtick', 2.^(0:9));
ylim([0,1].*ylim)
xlabel('Days since implant')
if P_in.ylabel_on
    switch P_in.metric
        case 'n_units'
            ylabel('Units')
        case 'n_good_units'
            ylabel('Single units')
    end
end
title(title_text, 'fontweight', 'normal')
text(1,30,['N = ' num2str(numel(rat_name)) ' rats'])