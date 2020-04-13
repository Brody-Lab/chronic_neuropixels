% PLOT_CHOICE_MOD Make a heat map showing the choice modulation
%
%=INPUT
%
%   time_s
%       An array specifying the time points
%
%   auROC
%       Area under an Receiver Operating Characteristic (ROC) classifying
%       between left and right choice trials. The values range from 0 to 1.
%       Values larger than 0.5 indicating a larger firing rates for right
%       choice trials
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   ax
%       The handle to an Axes object. If empty (default), then a figure is
%       created.
%
%   axes_properties
%       A cell array specifying the Axes properties
%
%   cell_area
%       A cell or string array specifying the location of the cells
%
%   color_bar
%       A scalar logical specifying whether to show a color bar
%
%   color_map_range
%       The data range over which the color indices are computed.
%
%   n_colors
%       Number of shades of colors
%
%=OPTIONAL OUTPUT
%
%   1) the handle to the colorbar
function [varargout] = plot_choice_mod(time_s, auROC, varargin)
P = inputParser;
addParameter(P, 'ax', [],  @(x) isa(x, 'matlab.graphics.axis.Axes') || isempty(x))
addParameter(P, 'axes_properties', {}, @(x) iscell(x))
addParameter(P, 'cell_area', {}, @(x) iscell(x)||isstring(x))
addParameter(P, 'color_bar', true,  @(x) isscalar(x) && islogical(x))
addParameter(P, 'color_map_range', [],  @(x) isnumeric(x)&&(isempty(x)||numel(x)==2))
addParameter(P, 'merge_LR', false, @(x) isscalar(x) && islogical(x))
addParameter(P, 'n_colors', 51, @(x) isscalar(x) && isnumeric(x))
addParameter(P, 'reference_event', {}, @(x) iscell(x)||ischar(x)||iscell(x))
parse(P, varargin{:});
P = P.Results;
P.color_map_range = sort(P.color_map_range);
if ~isempty(P.ax)
    axes(P.ax)
else
    tall_figure
end
fig_prepare_axes        
if P.merge_LR
    auROC = abs(auROC-0.5)+0.5;
    P.color_map_range = [0.5,1];
    cmap = colormap('gray');
    colormap(flip(cmap));
    color_index_center = 0.75;
else
    cmap = PB_get_my_colors(P.n_colors);
    colormap(cmap)
    color_index_center = 0.5;
end
n_colors = size(cmap,1);
[color_index, color_bin_centers] = get_color_index(auROC, n_colors, ...
                                                    'center', color_index_center, ...
                                                    'range', P.color_map_range);
im_hdl = image(time_s, ...
               1:size(auROC,1), ...
               color_index, 'CDataMapping','direct');
           
set(gca, 'YLim', [0, size(auROC,1)+1], ...
         'xlim', [min(time_s), max(time_s)], ...
         'yTick', [1, size(auROC,1)])
kText = PB_get_constant('text');
if ~isempty(P.reference_event)
    xlabel(['Time from ' kText.(P.reference_event) ' (s)'])
else
    xlabel('Time (s)')
end
if P.color_bar
    hdl = colorbar;
    hdl.Label.String = 'Choice selectivity';
    hdl.Limits = [1, n_colors];
    ticks = [1, round(n_colors/2), n_colors];
    hdl.Ticks = ticks;
    ticklabels = color_bin_centers(ticks);
    ticklabels = cellfun(@(x) num2str(x, '%0.2f'), num2cell(ticklabels), 'uni', 0);
    hdl.TickLabels = ticklabels;
    hdl.Location = 'westoutside';
end
% label brain areas
if ~isempty(P.cell_area)
    kColor = PB_get_constant('color');
    k = 0;
    for ca = unique(P.cell_area)'
        vl_cells = find(P.cell_area == ca);
        y_loc = (min(vl_cells) + max(vl_cells))/2;
        x_loc = min(xlim)-range(xlim)/2;
        k = k + 1;
        text(x_loc, y_loc, ca, 'Color', kColor.default(k,:));
        plot(xlim, min(vl_cells)*[1,1]-0.3, 'Color', kColor.default(k,:));
        plot(xlim, max(vl_cells)*[1,1]+0.3, 'Color', kColor.default(k,:));
    end
end

if nargout > 0
    if exist('hdl', 'var')
        varargout{1} = hdl;
    else
        varargout{1} = [];
    end
end