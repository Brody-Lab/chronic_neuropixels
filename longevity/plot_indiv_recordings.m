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
function [] = plot_indiv_recordings(Cells, brain_region, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'metric', 'n_units', ...
    @(x) validateattributes(x, {'string', 'cell', 'char'}, {}))
addParameter(parseobj, 'varying', 'marker', @(x) ismember(x, {'color', 'marker'}));
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
assert(ismember(brain_region, {'mFC', 'MCtx_ADS', 'AVS'}))
switch brain_region
    case 'mFC'
        rat_name = ["T212"; "T224"; "T176"];
        unique_bank = 0;
        title_text = 'Medial frontal ctx';
    case 'MCtx_ADS'
        rat_name = ["T181"; "T182"; "T219"];
        unique_bank = 1;
        title_text = 'Motor ctx and dorsal striatum';
    case 'AVS'
        rat_name = ["T181"; "T182"; "T219"];
        unique_bank = 0;
        title_text = 'Ventral striatum';
end
P = get_parameters;
T = make_table_from_Cells(Cells);
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
for i = 1:numel(rat_name)
    idx=T.rat == rat_name(i) & T.unique_bank == unique_bank;
    % if there are multiple probes per animal, use only the one with the
    % largest number of recordings
    Ct = groupcounts(T(idx,:), 'probe_serial'); 
    idx = idx & T.probe_serial == Ct.probe_serial(end);
    switch P_in.varying
        case 'color'
            hdl(i)= plot(T.days_since_surgery(idx)+1, T.(P_in.metric)(idx), ...
                 'o-', 'Color', P.color_order(i,:));
        case 'marker'
            hdl(i) = plot(T.days_since_surgery(idx)+1, T.(P_in.metric)(idx), ...
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
title(title_text)
text(1,30,['N = ' num2str(numel(rat_name)) ' rats'])