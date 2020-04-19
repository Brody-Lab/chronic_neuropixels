% FIT_EXP_DECAY_TO_DATA fit a sum of exponential decays to the data
%
%=INPUT
%
%   T
%       A table made from GET_METRICS_FROM_CELLS
%
%=OPTIONAL INPUT
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
%   normaliz_initial_value
%       Normalize the initial value
%       
function [p_hat, p_CI, p_boot] = fit_exp_decay_to_data(T, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))

addParameter(parseobj, 'n_boot',P.longevity_n_boots, ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 1}))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'normalize_initial_value', true, @(x) islogical(x)&&isscalar(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
%% normalize by electrodes
if numel(unique(T.condition)) > 1 && ~P_in.normalize_by_electrodes
    P_in.normalize_by_electrodes=true;
end
if P_in.normalize_by_electrodes && ~strcmp(P_in.metric, 'Vpp')
    metric=T.(P_in.metric)./T.n_elec;
else
    metric=T.(P_in.metric);
end
for i = 1:numel(unique(T.condition))
    [x,sort_idx] = sort(T.days_elapsed(T.condition==i));
    y = metric(T.condition==i);
    y=y(sort_idx);
    x=double(x);
    y=double(y);
    idx_nan = isnan(x)|isnan(y);
    x =x(~idx_nan);
    y=y(~idx_nan);
    if P_in.normalize_initial_value
        [p_hat{i} p_CI{i}, ~, p_boot{i}] = fit_sum_2_exp_decay_norm(x,y, 'n_boot', P_in.n_boot);
    else
        [p_hat{i} p_CI{i}, ~, p_boot{i}] = fit_sum_2_exp_decay(x,y, 'n_boot', P_in.n_boot);
    end
end