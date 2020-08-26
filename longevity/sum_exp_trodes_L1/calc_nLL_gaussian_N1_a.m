function nLL = calc_nLL_gaussian_N1_a(b, XN1, Xk, t, y)
% CALC_NLL_GAUSSIAN calculate the negative loglikelihood of the
% sum-of-exponental model assuming gaussian noise using the terms N1 and
% a. The constant term and a factor of 1/2 are ignored.
%
%=INPUT
%
%   b
%       A vector of betas
%
%   XN1f
%       The design matrix in the term N1f
%
%   XN1s
%       The design matrix in the term N1s
%
%   Xk
%       The design matrix in the term k
%
%=OUTPUT
%      
%   nLL
%       negative log-likelihood

    yhat = calc_resp_var_N1_a(b, XN1, Xk, t);
    nLL = sum((y-yhat).^2); % ignoring the term independent of lambda
end