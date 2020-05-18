% SUM_2_EXP_DECAY The sum of two exponential decays
%
%   y = p(1) * (p(2)*exp(x*p(3)) + (1-p(2))*exp(x*p(4)));
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
%
%=OPTIONAL INPUT
%
%   x0
%       The initial days, by default 1 and specified in GET_PARAMETERS
function y = sum_2_exp_decay(x,p, varargin)
P=get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
x0=P_in.x0;
y = p(1) * (p(2)*exp((x-x0)*p(3)) + (1-p(2))*exp((x-x0)*p(4)));