function nLL = calc_nLL_poisson_N1f_N1s(b, XN1f, XN1s, Xk, t, y)
% CALC_NLL_POISSON_N1F_N1S calculate the negative loglikelihood of the
% sum-of-exponental model assuming Poisson noise using the terms N1f and
% N1s. Constants are excluded from the design matrices.
%
%=INPUT
%
%   b
%       A vector of betas
%
%   XN1f
%       The design matrix in the term N1f, excluding the constant
%
%   XN1s
%       The design matrix in the term N1s, excluding the constant
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

    yhat = calc_resp_var_N1f_N1s(b, XN1f, XN1s, Xk, t);
    nLL = sum(yhat) - y'*log(yhat);
end