function nLL = calc_nLL_poisson_N1_a(b, XN1, Xk, t, y)
% CALC_NLL_POISSON_N1_a calculate the negative loglikelihood of the
% sum-of-exponental model assuming Poisson noise using the terms N1 and
% a. Constants are excluded from the design matrices.
%
%=INPUT
%
%   b
%       A vector of betas
%
%   XN1
%       The design matrix in the term N1, excluding the constant
%
%   Xk
%       The design matrix in the term k, excluding the two constants
%
%   t
%       Number of days since implant minus one (t=num_days-1)
%
%   y
%       A vector of the observed values of the response variable.
%=OUTPUT
%      
%   nLL
%       negative log-likelihood. The constant terms that are independent of
%       the parameters are not included in the sum. 

    yhat = calc_resp_var_N1_a(b, XN1, Xk, t);
    nLL = sum(yhat) - y'*log(yhat);
end