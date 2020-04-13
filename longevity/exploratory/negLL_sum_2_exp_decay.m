% negLL_SUM_2_EXP_DECAY Calculate the negative log-likelihood of the sum of
% two exponential decays given 1-D dependent and independent data and
% parameters, assuming constant Gaussian noise
%
% The model is 
%
%   Y = A*f*exp(-X/t1) + A*(1-f)*exp(-X/t2) + e
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
%       parameters, [A1, t1, A2, t2];
%
%=OUTPUT
%   LL
%       the loglikelihood
function negLL = negLL_sum_2_exp_decay(x,y,p)
y_hat = sum_2_exp_decay(x,p);
n = numel(y);
residuals = y-y_hat;
sigma2_hat = 1/(n-2) * sum(residuals.^2); % assume constant noise
LL = -n/2*log(2*pi) - n*log(sigma2_hat) - 1/(2*sigma2_hat)*sum(residuals.^2);
negLL = -LL;