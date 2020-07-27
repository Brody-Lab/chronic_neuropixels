% FIT_SUM_2_EXP_DECAY fit a sum of two exponentia decays using least
% squares with bounds on the parameters.
%
%   y = p(1)*p(2)*exp((x-x0)*p(3)) + p(1)*(1-p(2))*exp((x-x0)*p(4))
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
%   fit_initial_value
%       A scalar logical specifiying whether the initial value of the
%       response variable (i.e., p(1)) will be fit or assumed to be the
%       mean of the response variables at x==x0.
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
%   x0
%       The initial days, by default 1 and specified in GET_PARAMETERS
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
function [p_hat, p_CI, p_boot] = fit_sum_2_exp_decay(x,y, varargin)
    P=get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'fit_initial_value', true, ...
        @(x) isscalar(x)&&(x==0||x==1))
    addParameter(parseobj, 'lower_bounds', [min(y), 0, -1e6, -1e6], ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
    addParameter(parseobj, 'n_boot',P.exp_decay_n_boots, ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 1}))
    addParameter(parseobj, 'p0', [1,0.5,-0.5,-0.01], ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
    addParameter(parseobj, 'upper_bounds', [2*max(y), 1, 10, 10], ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 4}))
    addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), ...
        @(x) validateattributes(x, {'numeric'}, {'scalar'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    x=x(:);
    y=y(:);
    assert(numel(x)==numel(y), 'The number of elements of x and y must match')
    %% Optimization settings
    opts=optimoptions('fmincon');
    opts.Display = 'notify';
    opts.MaxFunctionEvaluations = 1e4;
    opts.MaxIterations = 1e4;
    %% remove days earlier than the specified initial day
    idx = x>=P_in.x0;
    x=x(idx);
    y=y(idx);
    if sum(x==P_in.x0)<1 && ~P_in.fit_initial_value
        error('There is no data corresponding to %i days after the surgery', P_in.x0);
    end
    %% initial value
    if ismember('p0',parseobj.UsingDefaults) && ...
       sum(x==P_in.x0)>0
        P_in.p0(1) = mean(y(x==P_in.x0));
    end
    if P_in.fit_initial_value
        y0 = [];
    else
        y0 = mean(y(x==P_in.x0));
    end
    %% Fit
    if P_in.n_boot > 0
        fprintf('\nFitting %i bootstrap draws...', P_in.n_boot); tic
            n = numel(x);
            ind_init = find(x==P_in.x0);
            ind_subseq = find(x>P_in.x0);
            n_init = numel(ind_init);
            n_subseq = numel(ind_subseq);
            p_boot = cell(P_in.n_boot, 1);
        parfor i =1:P_in.n_boot
            if P_in.fit_initial_value
                idx = datasample(1:n, n);
            else
                idx = [datasample(ind_init,n_init);datasample(ind_subseq,n_subseq)];
            end
            xboot = x(idx);
            yboot = y(idx);
            p_boot{i} = nan(1,4);
            if P_in.fit_initial_value
                [p_minNegLL,~,exit_flag] = ...
                    fmincon(@(p) f_negLL_4(xboot,yboot,p,P_in.x0,y0), ...
                                         P_in.p0, ...
                                        [0, 0,1,-1], 0, ... % t2-t1 >= 0, or t2>=t1
                                        [], [], ...
                                        P_in.lower_bounds, ...
                                        P_in.upper_bounds, ...
                                        [], ...
                                        opts);
                if exit_flag > 0
                    p_boot{i} = p_minNegLL;
                end
            else
                [p_minNegLL,~,exit_flag] = ...
                fmincon(@(p) f_negLL_3(xboot,yboot,p,P_in.x0,y0), ...
                                     P_in.p0(2:end), ...
                                    [0,1,-1], 0, ... % t2-t1 >= 0, or t2>=t1
                                    [], [], ...
                                    P_in.lower_bounds(2:end), ...
                                    P_in.upper_bounds(2:end), ...
                                    [], ...
                                    opts);
                if exit_flag > 0
                    p_boot{i}(2:4) = p_minNegLL;
                end
            end
        end
        p_boot = cell2mat(p_boot);
        p_boot = p_boot(all(~isnan(p_boot),2),:); % remove iterations that failed to fit
        fprintf(' took %0.f seconds\n', toc);
        p_hat = median(p_boot);
        p_CI = [quantile(p_boot, 0.025); quantile(p_boot, 0.975)];
    else
        [p_hat,~,~,~,~,~,hessian] = ...
            fmincon(@(p) f_negLL_4(x,y,p,P_in.x0,y0), ...
                    P_in.p0, ...
                    [0, 0,1,-1], 0, ... 
                    [], [], ...
                    P_in.lower_bounds, ...
                    P_in.upper_bounds, ...
                    [], ...
                    opts);
        if nargout> 1
            se = sqrt(diag(hessian^-1));
            se=se(:)';
            v95 = norminv(0.975);
            p_CI = p_hat + v95*[-1;1].*se;
            p_boot = [];
        end
    end
    if ~P_in.fit_initial_value
        p_hat(1) = y0;
        if nargout> 1
            p_CI(:,1) = y0*[1;1];
        end
        p_boot(:,1) = y0*ones(size(p_boot,1),1);
    end
end
%% f_negLL_4
% Calculate the negative log-likelihood of the sum of
% two exponential decays given 1-D dependent and independent data and
% parameters, assuming constant Gaussian noise
%
% The model is 
%
%   y = p(1)*p(2)*exp(x*p(3)) + p(1)*(1-p(2))*exp(x*p(4)) + e
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
%   x0
%       The initial day
%
%   y0
%       Initial values of the response variable will be fit. If not
%       supplied or empty, the initial value will be fit. 
%
%=OUTPUT
%   LL
%       the loglikelihood
function negLL = f_negLL_4(x,y,p,x0,y0)
    y_hat = (p(2)*exp((x-x0)*p(3)) + (1-p(2))*exp((x-x0)*p(4)));
    if nargin < 4 ||isempty(y0)
        y_hat = y_hat * p(1);
    else
        y_hat = y_hat * y0;
    end
    n = numel(y);
    residuals = y-y_hat;
    sigma2_hat = 1/(n-2) * sum(residuals.^2); % though homoskedasticity most likely not true
    LL = -n/2*log(2*pi) - n*log(sigma2_hat) - 1/(2*sigma2_hat)*sum(residuals.^2);
    negLL = -LL;
end
%%
function negLL = f_negLL_3(x,y,p,x0,y0)
    y_hat = y0*(p(1)*exp((x-x0)*p(2)) + (1-p(1))*exp((x-x0)*p(3)));
    n = numel(y);
    residuals = y-y_hat;
    sigma2_hat = 1/(n-2) * sum(residuals.^2); % though homoskedasticity most likely not true
    LL = -n/2*log(2*pi) - n*log(sigma2_hat) - 1/(2*sigma2_hat)*sum(residuals.^2);
    negLL = -LL;
end