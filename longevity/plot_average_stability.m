% PLOT_AVERAGE_STABILITY plot the neuronal stability metrics averaged
% across recording sessions-conditions. Errorbars indicate the standard
% error of the mean.
%
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
%=OPTIONAL OUTPUT
%   
%   hdl
%       The handles of the shadecolorbar objects.
function [varargout] = plot_average_stability(T, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'color_order_offset', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}))
addParameter(parseobj, 'fit_type', 'power', @(x) ischar(x))
addParameter(parseobj, 'legend_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'normalize_initial_value', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'FaceAlpha', 0.2, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}))
addParameter(parseobj, 'print_sample_size', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x));
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
y_label = P.text.(P_in.metric);
%% normalize by electrodes
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
if P_in.normalize_initial_value
    y_label(1) = lower(y_label(1));
    y_label = ['Norm. ' y_label];
end
%% Plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
condition_names = get_condition_names(T);
for i = 1:numel(unique(T.condition))
    fprintf('\n"%s" Condition %i, %s', P_in.metric, i, condition_names{i});
    bootstat = nan(P.longevity_n_boots,max(T.days_bin));
    for j=1:max(T.days_bin)
        idx = T.condition==i & T.days_bin==j;
        if sum(idx) > 1
            bootstat(:,j) = bootstrp(P.longevity_n_boots,@mean,metric(idx(:)));
        elseif sum(idx)==1
            bootstat(:,j) = metric(idx(:));
        else
            continue
        end
        N(j) = sum(idx); % smaple size
    end
    if P_in.normalize_initial_value
        y0 = mean(metric(T.condition==i & T.days_elapsed == P_in.x0));
        bootstat = bootstat/y0;
    end
    hdl(i)=shadedErrorBar(P.longevity_time_bin_centers, bootstat,{@nanmean,@nanstd});
    hold on;
    hdl(i).mainLine.LineWidth=1;
    hdl(i).mainLine.LineStyle = '-';
    hdl(i).mainLine.Marker = 'o';
    
    if numel(unique(T.condition)) == 1
        the_color = [0,0,0];
    else
        clr_idx = i + P_in.color_order_offset;
        the_color = P.color_order(clr_idx,:);
    end
    hdl(i).mainLine.Color = the_color;
    if ~isempty(hdl(i).patch)
        hdl(i).patch.FaceAlpha=P_in.FaceAlpha;
        hdl(i).patch.FaceColor = the_color;
    end
    % display power law fit
    switch P_in.fit_type
        case 'exponential'
            idx = T.condition==i & T.days_elapsed>=P_in.x0;
            [x,sort_idx] = sort(T.days_elapsed(idx));
            y = metric(idx);
            y=y(sort_idx);
            x=double(x);
            y=double(y);
            idx_nan = isnan(x)|isnan(y);
            x =x(~idx_nan);
            y=y(~idx_nan);
            p_hat = fit_sum_2_exp_decay(x,y, 'n_boot', 0);
            y_hat = sum_2_exp_decay(x,p_hat);
            if P_in.normalize_initial_value
                y_hat = y_hat/mean((y_hat(x == P_in.x0)));
            end
            plot(x, y_hat, ':', 'Color', the_color, 'linewidth', 2)
    end
    % display sample sizes
    if P_in.print_sample_size
        str_N = sprintf('%0.f, ',N);
        fprintf('\n    N = (%s)',  str_N(1:end-2));
    end
    % display ratio of loss
    early = metric(T.condition == i & T.days_elapsed <=1);
    steady_state = metric(T.condition == i & T.days_elapsed > 30);
    early_boot = bootstrp(1e4, @(x) mean(x), early(:));
    steady_boot = bootstrp(1e4, @(x) mean(x), steady_state(:)); %
    steady2early = steady_boot./early_boot;
    fprintf('\n    Ratio of %s >14 days after implantation to its initial value: %0.2f [%0.2f, %0.2f]', ...
    P_in.metric, median(steady2early), quantile(steady2early, 0.025), quantile(steady2early, 0.975))
end
fprintf('\n')
if numel(unique(T.condition))>1 && P_in.legend_on
    % flipping so that for comparing among DV positions, the most
    % superficial bin is listed first.
    legend(flip([hdl.mainLine]),flip(condition_names), 'location', 'best');
end
ylim([0,1].*ylim)
xlabel('Days since implant')
if P_in.ylabel_on
    ylabel(y_label); 
end
if nargout > 0
    varargout{1} = hdl;
end