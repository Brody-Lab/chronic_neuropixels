%   Calculates the negative loglikelihood of the parameters of a psychometric function in the form
%   of a four-parameter logistic function in generating a set of observed behavioral responses. The
%   frequency poking right at each click difference is assumed to be a sample of an independent
%   (questionable assumption) binomial random variable. So if there were twenty unique click
%   differences, then there are twenty indepndent binomial random variables.
%
%   This function used for maximum likelihood estimation of the parameters of the psychometric
%   function.
%
%   params(1) =             gamma0, or asymptotic left poke frequency
%   params(1) + params(2) = gamma1, asymptotic right poke frequency
%   params(2) =             beta, or sensitivity
%   params(3) =             alpha, or bias
%
%   F(x) = gamma0 + gamma1/(1+exp(-sens*(x-bias)
%
%   N           number of trials for each click difference
%   K           number of pokedR for each click difference
%   clickDiff   click differences
%
%   2017-02-09      Created

function [neg_loglike] = negLogLike_logistic4(params, N, K, clickDiff)

gamma0  = params(1);
gamma1  = params(2);
sens    = params(3);
bias    = params(4);
    
nDiff = numel(clickDiff);

loglike_eachDiff = nan(1,nDiff);
for d = 1:nDiff
    n = N(d);
    k = K(d);
    x = clickDiff(d);    
    
    p = gamma0 + gamma1/(1+exp(-sens*(x-bias)));
    
    if isinf(nchoosek(n,k))
        % calculate this manually
        log_n_k_factorial = 0;
        for i = 2:(n-k)
            log_n_k_factorial = log_n_k_factorial + log(i);
        end
        log_n_factorial_truncated = 0;
        for i = n:-1:(k+1)
            log_n_factorial_truncated = log_n_factorial_truncated + log(i);
        end
        log_nchoosek = log_n_factorial_truncated - log_n_k_factorial;
    else
        log_nchoosek = log(nchoosek(n,k));
    end
    loglike_eachDiff(d) = log_nchoosek + k*log(p) + (n-k)*log(1-p);
%     if isnan(loglike_eachDiff(d))
%         fprintf(2, 'NaN negative loglikelihood values\n')
%         keyboard
%     end
end

% ignore the infinite values
loglike_eachDiff(isinf(loglike_eachDiff)) = NaN;

neg_loglike = -1 * nansum(loglike_eachDiff);