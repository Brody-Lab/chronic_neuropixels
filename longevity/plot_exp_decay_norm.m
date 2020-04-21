%=INPUT
%
%   T
%       A table made from GET_METRICS_FROM_CELLS
%
%=OPTIONAL INPUT
%
%   anatom_bin_edges
%       If CONDITION_ON were 'ML','AP', or 'DV', then this increasing
%       vector specifies the bin edges for binning anatomical coodinates
%
%   axes
%       An AXES objects where the plot will be made. If empty(default),
%       then a new figureis created.
%
%   FaceAlpha
%       A scalar value between 0 and 1 that specifies the opacity of the
%       shading that represent the standard error of mean. When a large
%       number of conditions are plotted, this value can be set to 0 for
%       better clarity.
%
%   fit_type
%       A char array specifying the model type to fit the metric as a
%       function of days elapsed.
%
%   legend_on
%       A scalar logical indicating whether to show the legend
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   normalize_by_electrodes
%       A scalar logical specifying whether to normalize neuronal metrics
%       by the number of electrodes that satisfy the same condition. 
%
%   print_sample_size
%       Logical scalar indiating whether to show the sample sizes in the
%       command line
%
%   ylabel_on
%       A scalar logical specifying whether to show the ylabel
%       
function [] = plot_exp_decay_norm(T, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'color_order_offset', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}))
addParameter(parseobj, 'component', 'both', @(x) ismember(x, {'both', ...
                                                           'fast', ...
                                                           'slow'}))
addParameter(parseobj, 'initial_day', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
addParameter(parseobj, 'legend_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'p_boot', {}, @(x) iscell(x))
addParameter(parseobj, 'plot_data', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'xlabel_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'show_equation', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
y_label = P.text.(P_in.metric);
y_label(1)=lower(y_label(1));
y_label = ['Norm. # of ' y_label];
if max(T.condition)>1
    P_in.normalize_by_electrodes = true;
end
if P_in.normalize_by_electrodes
    T.(P_in.metric) = T.(P_in.metric)./T.n_elec;
    y_label = [y_label ' / electrode'];
end
%%
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
% set(gca, 'pos', get(gca, 'pos').*[1,1,1,0.5])
if P_in.plot_data
    set(gca, 'YLim', [0, 1.5])
else
    set(gca, 'YLim', [0, 1])
end
condition_names = get_condition_names(T);
x_plot = 2.^(P.longevity_log2_time_bin_edges(1):0.1:P.longevity_log2_time_bin_edges(end));
x0 = min(P.longevity_time_bin_centers);
x_plot(x_plot<x0) = [];
for i = 1:max(T.condition)
    % figure out the color
    if numel(unique(T.condition)) == 1
        the_color = [0,0,0];
    else
        clr_idx = i + P_in.color_order_offset;
        the_color = P.color_order(clr_idx,:);
    end
    
    % plot fits
    x = T.days_elapsed(T.condition==i);
    y = T.(P_in.metric)(T.condition==i);
    y = y/mean(y(x==x0));
    if isempty(P_in.p_boot)
        p_hat = fit_sum_2_exp_decay_norm(x,y, 'n_boot', 0); 
    else
        p_hat = median(P_in.p_boot{i});
    end
    y_plot = calc_exp_decay(x_plot, x0, p_hat, P_in.component);
    hdl(i) = plot(x_plot,y_plot, ':', 'Color', the_color, 'linewidth', 2);
    
    % plot confidence interval of fits
    if ~isempty(P_in.p_boot)
        y_boot = calc_exp_decay(x_plot, x0, P_in.p_boot{i}, P_in.component);
        y_CI = quantile(y_boot, [0.025; 0.975]);
        shadeplot(x_plot, y_CI(1,:), y_CI(2,:), 'Color', the_color)        
    end
    
    % plot data
    if P_in.plot_data
        bootstat = nan(P.longevity_n_boots, max(T.days_bin));
        days_bin=T.days_bin(T.condition==i);
        x_plot = P.longevity_time_bin_centers;
        for j=1:max(T.days_bin)
            idx = days_bin==j;
            if sum(idx) > 1
                bootstat(:,j) = bootstrp(P.longevity_n_boots,@mean,y(idx));
            elseif sum(idx)==1
                bootstat(:,j) = y(idx);
                x_plot(j) = x(idx); % if there is a single sample, don't bin it.
            end
        end
        errorbar(x_plot, mean(bootstat), std(bootstat),  'o', ...
                 'Color', the_color, 'linewidth', 1);
    end
end
if P_in.xlabel_on
    xlabel('Days since implant')
end
if P_in.ylabel_on
    ylabel(y_label)
end
if P_in.show_equation
    switch P_in.component
        case 'slow'
            str = '$$ (1-\alpha) e^{-(t-1)/\tau_{s}} $$';   
        case 'fast'
            str = '$$ \alpha e^{-(t-1)/\tau_{f}}$$';   
        case 'both'
            str = '$$ \alpha e^{-(t-1)/\tau_{f}} + (1-\alpha) e^{-(t-1)/\tau_{s}} $$';   
    end
    txt_hdl = text(2, max(ylim)*1.1,str,'Interpreter','latex', 'fontsize', get(gca, 'FontSize')*1.1);
%     center_text(txt_hdl,1);
    txt_hdl.Position(1) = (max(xlim)/min(xlim)/txt_hdl.Extent(3)*txt_hdl.Extent(1))^0.5;
end


if numel(unique(T.condition))>1 && P_in.legend_on
    % flipping so that for comparing among DV positions, the most
    % superficial bin is listed first.
    legend(flip(hdl),flip(condition_names), 'location', 'best');
end
end
%% CALC_EXP_DECAY
% calculate exponential decay
%
%   y = p(1)*exp(-(x-x0)/p(2)) + (1-p(1))*exp(-(x-x0)/p(2))
%
%=INPUT 
%   x
%       time data
%
%   x0
%       time offset
%
%   p
%       parameters
%
%   component
%       A char array speciifying whether both or one of the two components
%       are plotted
%
%=OUTPUT
%   
%   y
%       dependent data
function y = calc_exp_decay(x, x0, p, component)
    a = p(:,1);
    tau1 = p(:,2);
    tau2 = p(:,3);
    switch component
        case 'both'
            y = a.*exp(-(x-x0)./tau1) + (1-a).*exp(-(x-x0)./tau2);
        case 'fast'
            y = a.*exp(-(x-x0)./tau1);
        case 'slow'
            y = (1-a).*exp(-(x-x0)./tau2);
    end
end