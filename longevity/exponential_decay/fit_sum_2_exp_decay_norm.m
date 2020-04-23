% FIT_SUM_2_EXP_DECAY fit a sum of two exponentia decays using least
% squares with bounds on the parameters.
%
%   y = p(1)*exp(-(x-x0)/p(2)) + (1-p(1))*exp(-(x-x0)/p(2))
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
function [p_hat, p_CI, LL, p_boot, exit_flag] = fit_sum_2_exp_decay_norm(x,y, varargin)
    parseobj = inputParser;
    addParameter(parseobj, 'lower_bounds', [min(y), 0, 1e-6, 1e-6], ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
    addParameter(parseobj, 'p0', [1,1,1,100], ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
    addParameter(parseobj, 'upper_bounds', [2*max(y), 1, 1e6, 1e6], ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
    addParameter(parseobj, 'n_boot',0, ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 1}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    P=get_parameters;
    x=x(:);
    y=y(:);
    assert(numel(x)==numel(y), 'The number of elements of x and y must match')
    %% Optimization settings
    opts=optimset('fmincon');
    opts.Display = 'off';
    opts.TolFun = 1e-10;
    %% remove days earlier than the specified initial day
    x0= min(P.longevity_time_bin_centers);
    idx = x>=x0;
    x=x(idx);
    y=y(idx);
    %% normalize
    if sum(x==x0)<1
        error('There is data corresponding to %i days after the surgery', x0);
    end
    y_initial = mean(y(x==x0));
    y = y/y_initial;
    x = x - x0;
    p0 = [0.5, 0.5, 100];
    if P_in.n_boot > 0
        ind_init = find(x==x0);
        ind_subseq = find(x>x0);
        n_init = numel(ind_init);
        n_subseq = numel(ind_subseq);
        p_boot = nan(P_in.n_boot, 3);
        i = 0;
        fprintf('\nFitting %i bootstrap draws...', P_in.n_boot); tic
        while i <= P_in.n_boot
            idx = [datasample(ind_init,n_init);datasample(ind_subseq,n_subseq)];
            xboot = x(idx);
            yboot = y(idx);
            [p_minNegLL,~,exit_flag] = ...
                fmincon(@(p) f_negLL(xboot,yboot,p), ...
                                    p0, ...
                                    [0,1,-1], 0, ... % t1-t2 <= 0, or t1<=t2
                                    [], [], ...
                                    [0, 1e-6, 1e-6], ...
                                    [1, 1e6,  1e6], ...
                                    [], ...
                                    opts);
            if exit_flag > 0
                i = i+1;
                p_boot(i,:) = p_minNegLL;
            end
        end
        fprintf(' took %0.f seconds\n', toc);
        p_hat = median(p_boot);
        p_CI = [quantile(p_boot, 0.025); quantile(p_boot, 0.975)];
    else
        [p_hat,~,~,~,~,~,hessian] = ...
            fmincon(@(p) f_negLL(x,y,p), ...
                    p0, ...
                    [0,1,-1], 0, ... % t1-t2 <= 0, or t1<=t2
                    [], [], ...
                    [0, 1e-6, 1e-6], ...
                    [1, 1e6,  1e6], ...
                    [], ...
                    opts);
        se = sqrt(diag(hessian^-1));
        se=se(:)';
        v95 = norminv(0.975);
        p_CI = p_hat + v95*[-1;1].*se;
        p_boot = [];
    end
    LL = -f_negLL(x,y,p_hat);
end
%% f_negLL
% Calculate the negative log-likelihood of the sum of
% two exponential decays given 1-D dependent and independent data and
% parameters, assuming constant Gaussian noise
%
% The model is 
%
%   Y = f*exp(-X/t1) + (1-f)*exp(-X/t2) + e
%
% e ~ N(0,sigma2), independent of X, and Cov(ei,ej) = 0, i.e. independent
% across observations
%
%=INPUT
%   x
%       independent data, R^1
%
%   y
%       response data, R^1
%
%   p
%       parameters, [f, t1, t2]; "f" is the fraction of the population with
%       faster decline. "t1" is the time constant of the faster decay, and
%       "t2" is that of the slower decay.
%
%=OUTPUT
%   LL
%       the loglikelihood
function negLL = f_negLL(x,y,p)
    y_hat = (p(1)*exp(-x/p(2)) + (1-p(1))*exp(-x/p(3)));
    n = numel(y);
    residuals = y-y_hat;
    sigma2_hat = 1/(n-2) * sum(residuals.^2); % though homoskedasticity most likely not true
    LL = -n/2*log(2*pi) - n*log(sigma2_hat) - 1/(2*sigma2_hat)*sum(residuals.^2);
    negLL = -LL;
end