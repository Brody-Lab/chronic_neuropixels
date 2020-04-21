% FIT_EXP_DECAY Fit a data to a sum of exponential decays. The number of
% terms is determined by a likelihood ratio test
%
%=INPUT
%   x
%      independent data, 1-D
%
%   y
%       response data, 1-D
%
%=OPTIONAL INPUT
%
%   p_threshold
%       p-value threshold for statistical significance
%
%   p0
%       initial values
%
%=OUTPUT
%
%   p_hat
%       estimated values of the parameters 
%
%   p_CI
%       95% confidence intervals estimated using the hessian matrix
%
%   pval
%       p-value from the log-likelihood test
function [p_hat, p_CI] = fit_exp_decay(x,y, varargin)
parseobj = inputParser;
addParameter(parseobj, 'lower_bounds', [1, 0, 1e-6, 1e-6], ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
addParameter(parseobj, 'p_threshold', 0.01, ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 1, '>', 0, '<', 1}))
addParameter(parseobj, 'p0', [1,1,1,100], ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
addParameter(parseobj, 'upper_bounds', [inf, 1, 1e6, 1e6], ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
parse(parseobj, varargin{:});
P = parseobj.Results;
y0 = mean(y(x==min(x)));
[p_hat, p_CI] = fit_sum_2_exp_decay(x,y);
% fit the sum of an exponential and a constant
[p_hat{2}, p_CI{2}, LL(2)] = fit_sum_2_exp_decay(x,y, ...
                                        'p0', [y0,1,y0,inf], ...
                                        'lower_bounds', [0,0,0, 1e10], ...
                                        'upper_bounds', [y0_max,inf,y0_max, 1e10]);
% fit the sum of two exponential
[p_hat{3}, p_CI{3}, LL(3)] = fit_sum_2_exp_decay(x,y, ...
                                        'p0', [y0,1,y0,100], ...
                                        'lower_bounds', zeros(1,4), ...
                                        'upper_bounds', [y0_max,inf,y0_max,inf], ...
                                        'n_boot', 1000);                                    
AIC(1) = 2*(2-LL(1));
AIC(2) = 2*(3-LL(2));
AIC(3) = 2*(4-LL(3));
[~,idx] = min(AIC);
p_hat = p_hat{idx};
p_CI = p_CI{idx};