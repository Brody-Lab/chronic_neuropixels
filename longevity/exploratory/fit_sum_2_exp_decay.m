% FIT_SUM_2_EXP_DECAY fit a sum of two exponentia decays using least
% squares with bounds on the parameters.
%
%   y = p(1)*exp(-x/p(2)) + p(3)*exp(-x/p(4))
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
%   lower_bounds
%       lower bounds of the parameters during optimization
%
%   n_boot
%       Number of times for bootstrapping
%
%   p0
%       initial values
%
%   upper_bounds
%       upper bounds of the parameters during optimization
%
%=OUTPUT
%
%   p_hat
%       estimated values of the parameters 
%
%   p_CI
%       95% confidence intervals estimated using the hessian matrix
%
%   LL
%       The log-likelihood, for doing AIC, BIC, or log-likelihood test
%
%   hessian
%       The hessian matrix
function [p_hat, p_CI, LL, hessian] = fit_sum_2_exp_decay(x,y, varargin)
parseobj = inputParser;
addParameter(parseobj, 'lower_bounds', [1, 0, 1e-6, 1e-6], ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
addParameter(parseobj, 'p0', [1,1,1,100], ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
addParameter(parseobj, 'upper_bounds', [inf, 1, 1e6, 1e6], ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
addParameter(parseobj, 'n_boot',0, ...
    @(x) validateattributes(x, {'numeric'}, {'numel', 1}))
parse(parseobj, varargin{:});
P = parseobj.Results;
x=x(:);
y=y(:);
assert(numel(x)==numel(y), 'The number of elements of x and y must match')
opts=optimset('fmincon');
opts.Display = 'off';
opts.TolFun = 1e-10;
if P.n_boot > 0
    ind_early = find(x < 8);
    ind_steady=find(x>=8);
    n_early = numel(ind_early);
    n_steady = numel(ind_steady);
    p_boot = nan(P.n_boot, 4);
    for i = 1:P.n_boot
        fprintf('%i', i)        
        idx = [datasample(ind_early,n_early);datasample(ind_steady,n_steady)];
        xboot = x(idx);
        yboot = y(idx);
        y0 = mean(yboot(xboot==min(xboot)));
        p0 = [y0, 0.5, 0.5, 100];
        A_max = 2*max(yboot);
        A_min = 0.5*max(yboot);
        p_boot(i,:) = fmincon(@(p) negLL_sum_2_exp_decay(xboot,yboot,p), ...
                            p0, ...
                            [0,0,1,-1], 0, ... % t1-t2 <= 0, or t1<=t2
                            [], [], ...
                            [A_min, 0, 1e-6, 1e-6], [A_max, 1, 1e6, 1e6], ...
                            [], ...
                            opts);
    end
    
    p_hat = median(p_boot);
    p_CI = [quantile(p_boot, 0.025); quantile(p_boot, 0.975)];
else
    [p_hat,~,~,~,~,~,hessian] = ...
        fmincon(@(p) negLL_sum_2_exp_decay(x,y,p), ...
                P.p0, ...
                [], [], [], [], ...
                P.lower_bounds, P.upper_bounds, ...
                [], ...
                opts);
    se = sqrt(diag(hessian^-1));
    se=se(:)';
    v95 = norminv(0.975);
    p_CI = p_hat + v95*[-1;1].*se;
end
LL = -negLL_sum_2_exp_decay(x,y,p_hat);