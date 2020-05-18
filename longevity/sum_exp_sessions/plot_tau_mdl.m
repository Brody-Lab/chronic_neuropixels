function []=plot_tau_mdl(T_ed, varargin)
    P=get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
    addParameter(parseobj, 'color_order_offset', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}))
    addParameter(parseobj, 'condition', '', ...
        @(x) validateattributes(x, {'char'}, {}))
    addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
    addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && (x==0||x==1));      
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    
    idx = cellfun(@(x) contains(x, P_in.condition), T_ed.cond_name) & ...
          strcmp(T_ed.metric, P_in.metric);
    
    medians = T_ed.tau_mdl(idx);
    upper_err = T_ed.tau_mdl_CI(idx,2) - medians;
    lower_err = medians-T_ed.tau_mdl_CI(idx,1);
    if isempty(P_in.axes)
        figure('Position', P.figure_position_longevity)
    else
        axes(P_in.axes)
    end
    set(gca, P.axes_properties{:}, ...
             'XLim', [0.5, numel(medians)+0.5], ...
             'XTick', [], ...
             'YScale', 'log')
    for i = 1:numel(medians)
        errorbar(i, medians(i), lower_err(i), upper_err(i), 'o-', ...
                'color', P.color_order(i+P_in.color_order_offset,:), ...
                'linewidth', 1)
    end
    xlabel(P.text.(P_in.metric), 'fontweight', 'normal')
    if P_in. ylabel_on
        ylabel('Model time constant (days)')
    end
    ylim([10^0, 10^6])