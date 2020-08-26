function nLL = calc_nLL_poisson(b, XN1f, XN1s, Xk, t, y)
% CALC_NLL_POISSON calculate the negative loglikelihood of the
% sum-of-exponental model assuming Poisson noise. Constant terms are
% excluded
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
%     if any(yhat<0)
%         error('negative yhat?')
%     end
    nLL = sum(yhat) - y'*log(yhat); % ignoring the term independent of lambda
end