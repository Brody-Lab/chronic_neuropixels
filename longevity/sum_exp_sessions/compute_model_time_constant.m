% COMPUTE_MODEL_TIME_CONSTANT compute the effective time constant of the
% sum-of-exponentials model
%
%   y = p(1)*p(2)*exp(x*p(3)) + p(1)*(1-p(2))*exp(x*p(4))
%
%=INPUT
%
%   p
%       An N x 4 array. Each column of P corresponds to a parameter of the
%       model, and each row corresponds to an instance of the model
%
%=OUTPUT
%
%   t_const
%       A column vector with the same number of rows as p
%
%=OPTIONAL INPUT
%
%   x0
%       The initial days, by default 1 and specified in GET_PARAMETERS

function t_const = compute_model_time_constant(p, varargin)
P=get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
fprintf('\nComputing model time constant... '); tic
assert(size(p,2)==4, 'The array "p" needs to have four columns')
t_const = nan(size(p,1),1);
x = sym('x');
parfor i = 1:size(p,1)
    alpha = p(i,2);
    kf = p(i,3);
    ks = p(i,4);
    try
        eqn = exp(-1) == alpha*exp(kf*x) + (1-alpha)*exp(ks*x);
        t_const(i,1) = vpasolve(eqn, x);
    catch ME
        if ks > 0
            t_const(i,1) = 10^6;
        else
            rethrow(ME)
        end
    end
end
t_const = t_const + P_in.x0;
fprintf('took %0.f seconds\n', toc)