function nLL = calc_nLL_gaussian(b, XN1f, XN1s, Xk, t, y)
% CALC_NLL_GAUSSIAN calculate the negative loglikelihood of the
% sum-of-exponental model assuming gaussian noise. The oonstant term and a
% factor of 1/2 are ignored. 
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

    yhat = calc_resp_var(b, XN1f, XN1s, Xk, t);
    nLL = sum((y-yhat).^2); % ignoring the term independent of lambda
end