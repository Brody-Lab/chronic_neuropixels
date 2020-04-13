% SUM_2_EXP_DECAY The sum of two exponential decays
%
%   y = p(1) * (p(2)*exp(-x/p(3)) + (1-p(2))*exp(-x/p(4)));
%
%=INPUT
%   
%   x
%       A vector or a matrix
%
%   p
%       a four element array
%
%=OUTPUT
%
%   y
%       The response
function y = sum_2_exp_decay(x,p)
y = p(1) * (p(2)*exp(-x/p(3)) + (1-p(2))*exp(-x/p(4)));