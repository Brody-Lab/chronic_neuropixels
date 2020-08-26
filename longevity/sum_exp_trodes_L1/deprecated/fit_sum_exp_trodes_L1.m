% FIT_SUM_EXP_TRODES_L1 Fit the sum-of-exponentials model fit to the unit
% count of each electrode in each recording session, using a L1
% regularization.
%
% See
% ...\chronic_neuropixels\longevity\sum_exp_trodes_L1\model_specification.mlx
% for a description of the model.
%
%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILE and POSTPROCESS_CELLS
%
%=OUTPUT
%
%   S
%       A structure with the following fields
%       - T_mdl, a table indicating whether each parameter were fitted in
%       each model variant
%       - T_trode, a table of each electrode and each recoding
%       - T_regressor, a table specifying the range for each regressor
%       - T_res, a table of results
%       - P_in, a structure of the input parameters
%       - T_dsgn, the design matrix
%       
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   KFold
%       Number of cross-validation folds
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   noise
%       The noise model: "gaussian" or "poisson"
%
%   model_parameters
%       A char vector, string array, or a cell array of char that specified
%       the set of model parameters. Each parameter must be a member of
%       P.possible_model_parameters
%
%   num_phi
%       Number of phi values for each term (default: 25).
%
%   phi_ratio
%       Ratio of smallest to largest phi values (default: 1e-4). Phi is the
%       penalization paramter
%
%   shuffle
%       A logical scalar specifying whether to shuffle the rows for each
%       variable independently.
%       
function S = fit_sum_exp_trodes_L1(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.model_parameters_L1, ...
        @(x) all(ismember(x, P.sum_exp_trodes.model_parameters_L1)))
    addParameter(parseobj, 'noise', 'poisson', @(x) any(strcmpi(x, {'poisson', 'gaussian'})))
    addParameter(parseobj, 'num_phi', 12, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}))
    addParameter(parseobj, 'phi_ratio', 1e-4, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}))
    addParameter(parseobj, 'shuffle', false, @(x) x==0 || x==1)
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    [T_trode, T_regressor] = ...
              make_T_trode(Cells, 'model_parameters', P_in.model_parameters, ...
                                  'normalize_regressors', true, ...
                                  'unit_distance', P.unit_distance, ...
                                  'x0', P.x0);
    if P_in.shuffle
        n = size(T_trode,1);
        for i = 1:numel(T_regressor.name)
            T_trode.(T_regressor.name{i}) = T_trode.(T_regressor.name{i})(randperm(n));
        end
    end
    T_dsgn = make_design_matrix(T_trode, 'model_parameters', P_in.model_parameters);
    T_mdl = cell2table(num2cell(true(1,numel(P_in.model_parameters))));
    T_mdl.Properties.VariableNames = P_in.model_parameters;
    X = make_design_matrix_for_each_term(T_dsgn, T_mdl);
    S.X = X;
    
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
    n_N1f = size(X.N1f,2);
    n_N1s = size(X.N1s,2);
    n_k   = size(X.k,2);
    b0 = [-1; -0.01; ones(n_N1f,1); ones(n_N1s,1); -0.01*ones(n_k,1)]; 
       % [k_fast, k_slow, ...]
    [lb, ub] = get_bounds(X);
    [A,b] = get_linear_inequality_constraints(X);
          
    % randomly select half of the data for estimating phi
    n_trodes = size(X.N1f,1);
    i_est_phi = ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
    S.i_est_phi = i_est_phi;
    XN1f = X.N1f(i_est_phi,:);
    XN1s = X.N1s(i_est_phi,:);
    Xk = X.k(i_est_phi,:); 
    t = T_trode.days_since_init(i_est_phi);
    y = T_trode.(P_in.metric)(i_est_phi);
    
    %fit without penalization
    betas = fmincon(@(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y), b0, A, b, ...
                                      [], [], lb, ub, [], ...
                                      opts);
    B0.kf  = betas(1);
    B0.ks  = betas(2);
    B0.N1f = betas(3:n_N1f+2);
    B0.N1s = betas(n_N1f+3:n_N1f+n_N1s+2);
    B0.k   = betas(n_N1f+n_N1s+3:end);
    S.B0 = B0;
    
    % scaling parameter
    s = mean(abs([B0.N1f; B0.N1s]))/mean(abs(B0.k));
    beta_scale = [0; 0; ...
                  ones(n_N1f,1); ...
                  ones(n_N1s,1); ...
                  s*ones(n_k,1)];
    S.beta_scale = beta_scale;
    
    % Calculate the phi range. Range calculation is not yet principled: The
    % maximum phi value is exp(1) * sum(abs(B0.(term)));
    loghi = 0;
    loglo = log(P_in.phi_ratio);
    phi = exp(linspace(loghi,loglo,P_in.num_phi)+2);
    phi = abs(betas')*beta_scale*phi;
    n_mdl = size(phi,2);
    
    % fit with penalization
    cvp = cvpartition(sum(i_est_phi), 'KFold', P_in.KFold);
    betas_est_phi = cell(n_mdl,1);
    nLL_est_phi = cell(n_mdl, 1);
    tic
    parfor i = 1:n_mdl
        % cross-validate
        % To faciliate debugging, CROSSVAL is not used
        fprintf('\n%i', i)
        betas_est_phi{i} = nan(numel(b0), P_in.KFold);
        nLL_est_phi{i} = nan(1, P_in.KFold);
        f_c = @(x) abs(x')*beta_scale - phi(i);
        f_ceq = @(x) [];
        nonlinfcn = @(x) deal(f_c(x),f_ceq(x));
        for k = 1:P_in.KFold
            itrain = training(cvp,k);
            itest = test(cvp,k);    
            btrain = fmincon(@(betas) calc_nLL(betas, XN1f(itrain,:), XN1s(itrain,:), Xk(itrain,:), t(itrain), y(itrain)), ...
                             b0, [], [], [], [], [], [], ...
                             nonlinfcn, ...
                             opts);
            betas_est_phi{i}(:,k) = btrain;
            nLL_est_phi{i}(k) = calc_nLL(btrain, XN1f(itest,:), XN1s(itest,:), Xk(itest,:), t(itest), y(itest));
        end
    end
    toc
    S.betas_est_phi = betas_est_phi;
    S.nLL_est_phi = nLL_est_phi;
    
    % assemble a table including all the information for estimating phi
    T_phi.phi = phi(:);
    T_phi.LL_per_trode = cellfun(@(x) -1 * mean(x)/sum(i_est_phi), nLL_est_phi);
    for i = 1:n_mdl
        [~,idx] = min(abs(nLL_est_phi{i}-median(nLL_est_phi{i})));
        T_phi.betas(i,:) = betas_est_phi{i}(:, idx)';
    end
    T_phi = struct2table(T_phi);
    S.T_phi = T_phi;
    
    % estimate the parameter coefficients
    [~,i_opt_phi] = max(T_phi.LL_per_trode);
    S.i_opt_phi = i_opt_phi;
    XN1f = X.N1f(~i_est_phi,:);
    XN1s = X.N1s(~i_est_phi,:);
    Xk = X.k(~i_est_phi,:); 
    t = T_trode.days_since_init(~i_est_phi);
    y = T_trode.(P_in.metric)(~i_est_phi);
    f_c = @(x) abs(x')*beta_scale - phi(i_opt_phi);
    f_ceq = @(x) [];
    nonlinfcn = @(x) deal(f_c(x),f_ceq(x));    
    betas = fmincon(@(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y), ...
                             b0, [], [], [], [], [], [], ...
                             nonlinfcn, ...
                             opts);
    % collect variables into the output structure S
    S.P_in = P_in;
    S.T_trode = T_trode;
    S.T_regressor = T_regressor;
    S.T_dsgn = T_dsgn;
    S.T_mdl = T_mdl;
    S.betas = betas;
end
%% NONLCON
function [c,ceq] = nonlcon(betas, beta_scale, phi)
% NONLCON nonlinear constraints
%
%=INPUT
%
%   betas
%       A vector of model coefficients
%
%   beta_scale
%       A vector indicating the scaling of each beta
%
%   phi
%       A g-element vector indicating the penalization parameter of each
%       group
    c = abs(betas')*beta_scale - phi;
    ceq = [];
end 
%% GET_BOUNDS
%   Return lower and upper bounds of the parameter coefficients   
%
%=INPUT
%
%   X
%       The design matirx, a structure
%
%=OUTPUT
%
%   lb
%       Lower bounds. A real vector with the same number of elements as the
%       number of model parameters
%
%   ub
%       Upper bounds. A real vector with the same number of elements as the
%       number of model parameters
function [lb, ub] = get_bounds(X)
    n_N1f = size(X.N1f,2);
    n_N1s = size(X.N1s,2);
    n_k   = size(X.k,2);
    
    % lower and upper bounds
    lb = [-10; ... k_fast
          -1; ... k_slow
          -10*ones(n_N1f,1); ...
          -10*ones(n_N1s,1); ...
          -10*ones(n_k,1)];
    ub = [-1e-6; ... k_fast
          -1e-6; ... k_slow
          inf(n_N1f,1); ...
          inf(n_N1s,1); ...
          inf(n_k,1)];
end
%% GET_LINEAR_INEQUALITY_CONSTRAINTS
%   Return the linear inequality constraints expressed as an M-by-N matrix
%   A and an M-element vector b.
%
%=INPUT
%
%   X
%       The design matirx, a structure
%
%=OUTPUT
%
%   A
%       M by N matrix where M is the number of inequalities, and N is the
%       number of variables
%
%   b
%       M-element real vector
function [A,b] = get_linear_inequality_constraints(X)

    n_N1f = size(X.N1f,2);
    n_N1s = size(X.N1s,2);
    n_k   = size(X.k,2);

    % make k_fast <= k_slow
    A_kfs = [1,-1]; % k_fast - k_slow <= 0
    
    % make sure the N1s and N1f terms are positive for all regressors (who
    % are in the range of [0,1]). Negative because the constraint's format
    % is A*x <= b, and we want the linear combination to be >= 0
    if n_N1f > 1
        A_N1f = -double(dec2bin(0:2^(n_N1f-1)-1)=='1');
        A_N1f = [-ones(size(A_N1f,1),1), A_N1f]; % the constant term
    else
        A_N1f = -1; % the constant term
    end
    if n_N1s > 1
        A_N1s = -double(dec2bin(0:2^(n_N1s-1)-1)=='1');
        A_N1s = [-ones(size(A_N1s,1),1), A_N1s]; % the constant term
    else
        A_N1s = -1;
    end
    A = blkdiag(A_kfs, A_N1f, A_N1s);
    A = padarray(A, [0,n_k], 0, 'post');
    b = zeros(size(A,1),1);
end