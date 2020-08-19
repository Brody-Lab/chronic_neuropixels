% FIT_MDL_SUM_EXP_TRODES fit the sum-of-exponentials model to electrodes
%
%=INPUT
%
%   XN1
%       Design matrix, a structure
%
%   t
%
%       days from initial day
%
%   y
%
%       the response variable (unit count, etc.)
%=OPTIONAL INPUT
%
%   num_phi
%       Number of phi values for each term (default: 25).
%
%   phi_ratio
%       Ratio of smallest to largest phi values (default: 1e-4). Phi is the
%       penalization paramter
function [B, FitInfo] = fit_SoE_mdl_L1(X, t, y, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'noise', 'poisson', @(x) any(strcmpi(x, {'poisson', 'gaussian'})))
    addParameter(parseobj, 'num_phi', 25, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}))
    addParameter(parseobj, 'phi_ratio', 1e-4, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    FitInfo.P_in = P_in;
    switch P_in.noise
        case 'gaussian'
            calc_nLL = @calc_nLL_gaussian;
        case 'poisson'
            calc_nLL = @calc_nLL_poisson;
        otherwise
            error('unrecognized noise distribiution');
    end
    opts=optimoptions('fmincon');
    opts.Display = 'off';

    % Initial valuess
    n_XN1f = size(X.N1f,2);
    n_XN1s = size(X.N1s,2);
    n_Xk   = size(X.k,2);
    b0 = [-1; -0.01; ones(n_XN1f,1); ones(n_XN1s,1); -0.01*ones(n_Xk,1)]; 
       % [k_fast, k_slow, ...]
    
    % make k_fast <= k_slow
    A = zeros(size(b0));
    A = A(:)';
    A(2:3) = [1,-1]; % k_fast - k_slow <= 0
    b = 0;
        
    % lower and upper bounds
    lb = [-10; ... k_fast
          -1; ... k_slow
          -10*ones(n_XN1f,1); ...
          -10*ones(n_XN1s,1); ...
          -10*ones(n_Xk,1)];
    ub = [-1e-6; ... k_fast
          -1e-6; ... k_slow
          inf(n_XN1f,1); ...
          inf(n_XN1s,1); ...
          inf(n_Xk,1)];
    %fit without penalization
    betas = fmincon(@(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y), b0, A, b, ...
                                      [], [], lb, ub, [], ...
                                      opts);
    B0.kf  = betas(1);
    B0.ks  = betas(2);
    B0.N1f = betas(3:n_XN1f+2);
    B0.N1s = betas(n_XN1f+3:n_XN1f+n_XN1s+2);
    B0.k   = betas(n_XN1f+n_XN1s+3:end);
    FitInfo.B0 = B0;
    
    % not very principled: The maximum phi value is exp(1) * sum(abs(B0.(term)));
    loghi = log(exp(1));
    loglo = log(P_in.phi_ratio*exp(1));
    phi_range = exp(linspace(loghi,loglo,P_in.num_phi));
    for term = {'N1f', 'N1s', 'k'}
        term=term{:};
        Phi.(term) = phi_range * sum(abs(B0.(term)));
    end
end  
%% NONLCON
% 
%   nonlinear constraint
function [c,ceq] = nonlcon(betas, param_index, phi)

    for i = 1:numel(

end