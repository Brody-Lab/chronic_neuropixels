% STAT_TEST_STABILITY reports the p-value that a stability metric depends
% on the days elapsed for data in a specified day range
%
%=INPUT
%
%   T
%       A table made using GET_METRICS_FROM_CELLS
%
%=OUTPUT
%
%   pval
%       The probability (based on the t-statistic) that the slope of the fitted
%       line is zero. The i-th element corresponds to the i-th condition in
%       T
%
%=OPTIONAL INPUT
%
%   day_range
%       A two element vector that specifies the range of days elapsed
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform

function pval = stat_test_stability(T, varargin)
parseobj = inputParser;
addParameter(parseobj, 'day_range', [8, inf], @(x) validateattributes(x, {'numeric'}, {'numel', 2, 'increasing'}))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
idx_in_range = discretize(T.days_elapsed, P_in.day_range);
T = T(idx_in_range==1, :);
metric = T.(P_in.metric)./T.n_elec;
for i = 1:numel(unique(T.condition))
    idx = T.condition==i;
    pval(i,1)=anovan(metric(idx), T.days_elapsed(idx), 'display', 'off', 'continuous', 1);
end