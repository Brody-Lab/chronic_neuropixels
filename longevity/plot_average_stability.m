% PLOT_AVERAGE_STABILITY plot the neuronal stability metrics averaged
% across recording sessions-conditions. Errorbars indicate the standard
% error of the mean.
%
%=INPUT
%
%   Cells
%       A cel array of structures made by COLLECT_CELL_FILES and
%       POSTPROCESS_CELLS
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
%   condition_on
%       A char array specifying the different conditions used for
%       averaging. The option are {'bank', 'ML', 'AP', 'DV'}.
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
%       by the number of electrodes that satisfy the same condition. This
%       is should be turned on CONDITION_ON is not empty
%
%   print_sample_size
%       Logical scalar indiating whether to show the sample sizes in the
%       command line
%
%   ylabel_on
%       A scalar logical specifying whether to show the ylabel
%       
function [] = plot_average_stability(T, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'color_order_offset', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}))
addParameter(parseobj, 'fit_type', 'power', @(x) ischar(x))
addParameter(parseobj, 'legend_on', true, @(x) isscalar(x) && islogical(x));
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'FaceAlpha', 0.2, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}))
addParameter(parseobj, 'print_sample_size', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'ylabel_on', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;
t_bin_edges = P.longevity_time_bin_edges;
t_bin_centers = P.longevity_time_bin_centers;
n_boots = P.longevity_n_boots;
y_label = P.text.(P_in.metric);
% T=T(T.days_elapsed>0,:);
% T.days_elapsed = T.days_elapsed + 1;
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
%% Model
switch P_in.fit_type
    case 'power'
        modelfun = @(b,x) b(1)*x.^b(2);
        b0=[1,-0.5, 1, -0.1];
        if ~P_in.normalize_by_electrodes
            b0(1)=600;
        end
    case 'exponential'
        modelfun = @(b,x) b(1)*exp(b(2)*(x)) + b(3); % + b(3)*exp(b(4)*x);
        b0=[10,-2,0.5,-0.001];
    case 'log'
        modelfun = @(b,x) b(1) + b(2)*log(x) + b(3) + b(4)*log(x);
        b0=[1, -2, 1,-2];
end
if ~P_in.normalize_by_electrodes
    b0([1,3])=300;
end
opts = statset('nlinfit');
%% Plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
condition_names = get_condition_names(T);
clear boots
for i = 1:numel(unique(T.condition))
    fprintf('\n"%s" Condition %i, %s', P_in.metric, i, condition_names{i});
    boots = nan(n_boots,numel(t_bin_centers));
    for j=1:numel(t_bin_centers)
        idx = T.days_elapsed>=t_bin_edges(j) & T.days_elapsed<t_bin_edges(j+1) & ...
              T.condition==i;
        if sum(idx) > 1
            boots(:,j) = bootstrp(n_boots,@mean,metric(idx(:)));
        elseif sum(idx)==1
            boots(1:n_boots,j) = metric(idx(:));
        else
            continue
        end
        N(j) = sum(idx); % smaple size
    end
    hdl(i)=shadedErrorBar(t_bin_centers, boots,{@nanmean,@nanstd});
%     hdl(i)=shadedErrorBar(t_bin_centers,avg, err);
    hold on;
    hdl(i).mainLine.LineWidth=1;
    hdl(i).mainLine.LineStyle = '-';
    hdl(i).mainLine.Marker = 'o';
    hdl(i).patch.FaceAlpha=P_in.FaceAlpha;
    if numel(unique(T.condition)) == 1
        the_color = [0,0,0];
    else
        clr_idx = i + P_in.color_order_offset;
        the_color = P.color_order(clr_idx,:);
    end
    hdl(i).mainLine.Color = the_color;
    hdl(i).patch.FaceColor = the_color;
    % display power law fit
    switch P_in.fit_type
        case 'exponential'
            [x,sort_idx] = sort(T.days_elapsed(T.condition==i));
            y = metric(T.condition==i);
            y=y(sort_idx);
            x=double(x);
            y=double(y);
            idx_nan = isnan(x)|isnan(y);
            x =x(~idx_nan);
            y=y(~idx_nan);
            [p_hat, p_CI] = fit_exp_decay(x,y);
            plot(x, sum_2_exp_decay(x,p_hat), '--', 'Color', the_color, 'linewidth', 1)
            for j = 1:4
                fprintf('\n        p(%i) = %0.2f [%0.2f, %0.2f]', ...
                        j, p_hat(j), p_CI(1,j), p_CI(2,j))
            end
        case {'power', 'log'}
            x = T.days_elapsed(T.condition==i);
            y = metric(T.condition==i);
            beta = nlinfit(x,y, modelfun, b0, 'errormodel', 'proportional');
            x = 2.^(-1:0.1:10);
            plot(x, modelfun(beta, x), '--', 'Color', the_color, 'linewidth', 1)
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