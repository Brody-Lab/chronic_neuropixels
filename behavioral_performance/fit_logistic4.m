% FIT_LOGISTIC4 Maximum-likelihood estimates of the parameters of a
% logistic function
%
%   y = gamma0 + gamma1./(1 + exp(-sens*(x-bias)));
%
%   Implemented using FMINCON
%
%=INPUT
%       y
%           response data, vector
%
%       x
%           independent data, vector
%
%=OPTIONAL INPUT
%
%   lower_bound
%       A 1x4 vector specifying the lower bound for the parameters GAMMA0,
%       GAMMA1, SENS, and BIAS, respectively.
%
%   upper_bound
%       A 1x4 vector specifying the upper bound for the parameters GAMMA0,
%       GAMMA1, SENS, and BIAS, respectively.    
%
%   initial_values
%       A 1x4 vector specifying the initial values for the parameters
%       GAMMA0, GAMMA1, SENS, and BIAS, respectively.
function [Fit, varargout] = fit_logistic4(y, x, varargin)
    x = x(:);
    y = y(:);
    P = inputParser;
    addParameter(P, 'lower_bound', [0, 0, 0, min(x)], @isnumeric) % lapse [0,1], sensitivity [0 5], bias [min(X), max(X)]
    addParameter(P, 'upper_bound', [1, 1, 5, max(x)], @isnumeric)
    addParameter(P, 'initial_values',[0.1, 0.9, 0.5, 0], @(x) isnumeric(x))% a bit of lapse, mediocre sensitivity, no bias 
    parse(P, varargin{:});
    P = P.Results;
    X_unique = unique(x);
    n_X_unique = numel(X_unique);
    K = nan(n_X_unique,1);
    N = nan(n_X_unique,1);
    for c = 1:n_X_unique
        idx = X_unique(c) == x;
        K(c,1) = sum(y(idx)); 
        N(c,1) = sum(idx);
    end
    warning('off', 'MATLAB:nchoosek:LargeCoefficient')
    nLL0 = negLogLike_logistic4(P.initial_values, N, K, X_unique);
    if isinf(nLL0) || isnan(nLL0)
        error('Failure to evaluate initial values')
    end
    [x_fmin, f_fmin, exitflag, output, ~, grad, hessian] = ...
        fmincon(@(param) negLogLike_logistic4(param, N, K, X_unique), ...
                P.initial_values, ...
                [], [], [], [], ...
                P.lower_bound, P.upper_bound, ...
                @boundcon, ...
                optimset('DiffMinChange', 0.0001, ...
                         'MaxIter', 200, ...
                         'Algorithm', 'interior-point', ...
                         'Display', 'off')     );
    Fit = struct;
    Fit.gamma0 = x_fmin(1);
    Fit.gamma1 = x_fmin(2);
    Fit.sens   = x_fmin(3);
    Fit.bias   = x_fmin(4);

    est_var = diag(inv(hessian)); % estimated variance
    est_var(est_var < 0) = NaN; % hessian matrix not positive semidefinite
    est_std = est_var.^0.5;

    Fit.std.gamma0 = est_std(1);
    Fit.std.gamma1 = est_std(2);
    Fit.std.sens = est_std(3);
    Fit.std.bias = est_std(4);
    
    % 95% confidence interval
    tcv = tinv(0.975, numel(y)-2); % T critical value
    Fit.error95.gamma0 = Fit.std.gamma0 * tcv;
    Fit.error95.gamma1 = Fit.std.gamma1 * tcv;
    Fit.error95.sens = Fit.std.sens * tcv;
    Fit.error95.bias = Fit.std.bias * tcv;

    % *** vargout ***
    varargout{1} = f_fmin;
    varargout{2} = exitflag;
    varargout{3} = output;
    varargout{4} = grad;
    varargout{5} = hessian;
end

function [c,ceq] = boundcon(x)
    c(1) = x(1) + x(2) - 1; % gamma0 + gamma1  <= 1
    c(2) = -(x(1) + x(2));  % gamma0 + gamma1 >= 0
    ceq = [];
end