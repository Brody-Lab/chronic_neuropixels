function [] = make_Vpp_histogram(T, varargin)
% MAKE_VPP_HISTOGRAM make a histogram of the average peak-to-peak
% amplitudes
%
%=INPUT
%
%   T
%       A table made by "get_metrics_from_Cells"
%
%=OPTIONAL INPUT
%
%   bin_edges
%       The edges of the histogram bins
%
%   color_order
%       An n-by-3 array with values >= 0 and <= 1
%
%   color_order_offset
%       An integer that shifts the color order specifying the color
%       representing each group of sessions
%
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'bin_edges', 0:25:800, ...
        @(x) validateattributes(x, {'numeric'}, {'increasing', 'nonnegative'}))
    addParameter(parseobj, 'color_order', P.color_order, ...
        @(x) validateattributes(x, {'numeric'}, {'>=' 0, '<=', 1, 'ncols', 3}))
    addParameter(parseobj, 'color_order_offset', 0, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    
    for i = 1:numel(unique(T.condition))
        if numel(unique(T.condition)) == 1
            the_color = [0,0,0];
        else
            clr_idx = i + P_in.color_order_offset;
            the_color = P_in.color_order(clr_idx,:);
        end
        set(gca, P.axes_properties{:})
        h = histogram(vertcat(T.Vpp_distrib{T.condition==i}), P_in.bin_edges, 'normalization', 'probability');
        set(h, 'FaceColor', the_color, 'FaceAlpha', 0.2, 'EdgeColor', the_color, 'linewidth', 1)
        xlabel('Peak-to-peak amplitude (\muV)')
        ylabel('Fraction of single units')
        yticks(0:0.05:2);
        yticklabels(arrayfun(@num2str, yticks, 'uni',0))
    end
end