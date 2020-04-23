% COMPUTE_MODEL_TIME_CONSTANT compute the half-life of the a model of the sum
% of two exponential decay terms
%
%   y = p(1)*p(2)*exp(-x/p(3)) + p(1)*(1-p(2))*exp(-x/p(4))
%
%=INPUT
%
%   p
%       An N x 4 array. Each column of P corresponds to a parameter of the
%       model, and each row corresponds to an instance of the model
%
%=OUTPUT
%
%   t_half
%       A column vector with the same number of rows as p
%
%=OPTIONAL INPUT
%
%   x0
%       The initial days, by default 1 and specified in GET_PARAMETERS

function t_half = compute_model_time_constant(p, varargin)
P=get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
fprintf('\nComputing model half life... '); tic
assert(size(p,2)==4, 'The array "p" needs to have four columns')
syms x
t_half = nan(size(p,1),1);
for i = 1:size(p,1)
    alpha = p(i,2);
    tau_f = p(i,3);
    tau_s = p(i,4);
    eqn = exp(-1) == alpha*(exp(-1/tau_f))^x + (1-alpha)*(exp(-1/tau_s))^x;
    t_half(i,1) = vpasolve(eqn, x);
end
t_half = t_half + P_in.x0;
fprintf('took %0.f seconds\n', toc)