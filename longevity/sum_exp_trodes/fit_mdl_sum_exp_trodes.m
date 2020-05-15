% FIT_MDL_SUM_EXP_TRODES fit the sum-of-exponentials model to electrodes
%
%=INPUT
%
%   XN1
%       Design matrix in the N1 term
%
%   Xk
%
%       Design matrix in the decay rate term
%
%
%   t
%
%       days from initial day
%
%   y
%
%       the response variable (unit count, etc.)
function b = fit_mdl_sum_exp_trodes(XN1, Xk, t, y, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'noise', 'gaussian', @(x) any(strcmpi(x, {'poisson', 'gaussian'})))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    switch P_in.noise
        case 'gaussian'
            calc_err = @mean_square_error;
        case 'poisson'
            calc_err = @poisson_nLL_approx;
        otherwise
            error('unrecognized noise distribiution');
    end
    opts=optimoptions('fmincon');
    opts.Display = 'off';
%     opts.MaxFunctionEvaluations = 1e5;
%     opts.MaxIterations = 1e4;
    % Initial valuess
    n_XN1 = size(XN1,2);
    n_Xk = size(Xk,2);
    b0 = [0.5; -1; -0.01; ones(n_XN1,1); -0.01*ones(n_Xk,1)]; % [alpha, tau_alpha, ...]
    % The linear constraints matrix
    % A*b >= 0
    % The weighted sum of the decay terms has to be positive.
    % All values of X are normalized to be in [0,1]
    if n_Xk > 1
        A = double(dec2bin(0:2^(n_Xk-1)-1)=='1');
        A = [ones(size(A,1),1), A]; % the constant term
    else
        A = 1;
    end
    A = [A;-A]; % greater than, and then less than
    % add some columns of zeroes
    n_col_A = size(A,1);
    A = [zeros(n_col_A, 2), zeros(n_col_A, n_XN1), A];
    bound = [-1e-6*ones(n_col_A/2,1); ... A <= -10^-6 forces kappa to be negative
              1e6 *ones(n_col_A/2,1)]; %  A >= -10^6  but not too negative
    % combine the two linear inequality matrices
    lb = [0; ... alpha
          -10; ... k_fast
          -1; ... k_slow
          -10*ones(n_XN1+n_Xk,1)];
    ub = [1; ... alpha
          -1e-6; ... k_fast
          -1e-6; ... k_slow
          inf(n_XN1+n_Xk,1)];
    % make k_fast <= k_slow
    A = zeros(size(b0));
    A = A(:)';
    A(2:3) = [1,-1]; % k_fast - k_slow <= 0
    bound = 0;
    b = fmincon(@(b) calc_err(b, XN1, Xk, t, y), b0, A, bound, ...
                                [], [], lb, ub, [], ...
                                opts);
%     poisson_nLL_approx(b, XN1, Xk, t, y)
%     lambda = calc_lambda(b, XN1, Xk, t);
% %     close all
%     figure('Pos', [500, 400, 500,  420])
%     plot(lambda, 'o')
%     set(gca, 'ylim', [0, 2], 'xtick', [])
%     hold on
%     refline(0,1)
%     ylabel('\lambda')
%     figure('Pos', [1100, 400, 500,  420])
%     LL = -lambda + y.*log(lambda) -log(factorial(y));
%     plot(LL, 'o')
%     set(gca, 'ylim', [-100, 0], 'xtick', [])
%     ylabel('LL')
%     title(['mean LL = ' num2str(mean(LL))])
end
%% POISSON_LL_APPROX
% Calculating the sum of the poisson loglikelihood given the observations
% and the parameters for approximating the lambda. The term in the
% loglikelihood that depends only on the observations is not included.
function nLL = poisson_nLL_approx(b, XN1, Xk, t, y)
    yhat = predict_y(b, XN1, Xk, t);
    nLL = sum(yhat) - sum(y.*log(yhat)); % ignoring the term independent of lambda
end
%% sum of absolute deviation
function mse = mean_square_error(b, XN1, Xk, t, y)
    yhat = predict_y(b, XN1, Xk, t);
    mse = mean((y-yhat).^2);
end