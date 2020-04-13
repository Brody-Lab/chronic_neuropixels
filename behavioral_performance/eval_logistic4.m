% EVAL_LOGISTIC4 evaluate a four-parameter logistic function
%
% y = gamma0 + gamma1./(1 + exp(-sens*(x-bias)));
%
%=INPUT
%   
%   x
%       input data
%
%   gamma0
%       lapse rate for trials for which the correct choice is left is
%       GAMMA0. [0-100], in percent.
%
%   gamma1
%       lapse rate parameter. The lapse rate for trials for which the
%       correct side is right is 100-(GAMMA0+GAMMA1)
%
%   sens
%       The slope of the logistic function at its midpoint
%
%   bias
%       The x-value for which y = GAMMA0 + GAMMA1 / 2
%       If GAMMA0 == 0 and GAMMA1 == 100, then this is the x-value for
%       which y = 50
%
function [y] = eval_logistic4(x, gamma0, gamma1, sens, bias)
    y = gamma0 + gamma1./(1 + exp(-sens*(x-bias)));