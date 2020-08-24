% fit_SOE_L1 Fit the sum-of-exponentials model fit to the unit
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
%   generate_data
%       Make fake response data using the observed predictors
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
%   num_lambda
%       Number of lambda values for each term (default: 25).
%
%   lambda_ratio
%       Ratio of smallest to largest lambda values (default: 1e-4). Lambda is the
%       penalization paramter
%
%   shuffle
%       A logical scalar specifying whether to shuffle the rows for each
%       variable independently.
%       
function S = fit_SoE_L1(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'generate_data', false, @(x) x==0 || x==1)
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.model_parameters_L1, ...
        @(x) all(ismember(x, P.sum_exp_trodes.model_parameters_L1)))
    addParameter(parseobj, 'noise', 'poisson', @(x) any(strcmp(x, {'poisson', 'gaussian'})))
    addParameter(parseobj, 'num_lambda', 12, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}))
    addParameter(parseobj, 'lambda_ratio', 1e-2, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}))
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
    
    if P_in.generate_data
        [ygen, bgen] = generate_response_var(X,  T_trode.days_since_init);
        T_trode.(P_in.metric) = ygen;
        S.bgen = bgen;
    end
    
    switch P_in.noise
        case 'gaussian'
            calc_nLL = @calc_nLL_gaussian;
        case 'poisson'
            calc_nLL = @calc_nLL_poisson;
    end
    
    opts=optimoptions('fmincon');
    opts.Display = 'notify';
    opts.MaxFunctionEvaluations = 1e5;
    % Initial valuess
    n_N1f = size(X.N1f,2);
    n_N1s = size(X.N1s,2);
    n_k   = size(X.k,2);
    b0 = [-1; -0.01; ones(n_N1f,1); ones(n_N1s,1); -0.01*ones(n_k,1)]; 
       % [k_fast, k_slow, ...]
%     [lb, ub] = get_bounds(X);
    [A,b] = get_linear_inequality_constraints(X);
          
    % randomly select half of the data for estimating lambda
    n_trodes = size(X.N1f,1);
    i_est_lambda = ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
    S.i_est_lambda = i_est_lambda;
    XN1f = X.N1f(i_est_lambda,:);
    XN1s = X.N1s(i_est_lambda,:);
    Xk = X.k(i_est_lambda,:); 
    t = T_trode.days_since_init(i_est_lambda);
    y = T_trode.(P_in.metric)(i_est_lambda);
    
    %fit without penalization
    tic
    betas = fmincon(@(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y), b0, A, b, ...
                                      [], [], [], [], [], ...
                                      opts);
    toc
    
    B0.kf   = betas(1);
    B0.ks   = betas(2);
    B0.N1f0 = betas(3);
    B0.N1f  = betas(4:n_N1f+2);
    B0.N1s0 = betas(n_N1f+3);
    B0.N1s  = betas(n_N1f+4:n_N1f+n_N1s+2);
    B0.k    = betas(n_N1f+n_N1s+3:end);
    assert(numel(cell2mat(struct2cell(B0))) == numel(b0));
    S.B0 = B0;
    
    % scaling parameter
    s = mean(abs([B0.N1f; B0.N1s]))/mean(abs(B0.k));
    beta_scale = [0; 0; 0; 
                  ones(n_N1f-1,1); ...
                  0;
                  ones(n_N1s-1,1); ...
                  s*ones(n_k,1)];
    S.beta_scale = beta_scale;
    
    nLL0 = calc_nLL(betas, XN1f, XN1s, Xk, t, y);
    loghi = -4; % this is empirical
    loglo = loghi+log(P_in.lambda_ratio);
    lambda = exp(linspace(loghi,loglo,P_in.num_lambda-1))*nLL0/(abs(betas')*beta_scale);
    lambda = [lambda,0];
    n_mdl = size(lambda,2);
     
    % fit with penalization
    cvp = cvpartition(sum(i_est_lambda), 'KFold', P_in.KFold);
    betas_est_lambda = cell(n_mdl,1);
    nLL_est_lambda = cell(n_mdl, 1);
    tic
    parfor i = 1:n_mdl
        % cross-validate
        % To faciliate debugging, CROSSVAL is not used
        fprintf('\n%i', i)
        betas_est_lambda{i} = nan(numel(b0), P_in.KFold);
        nLL_est_lambda{i} = nan(1, P_in.KFold);
        for k = 1:P_in.KFold
            itrain = training(cvp,k);
            itest = test(cvp,k);    
            fun = @(betas) calc_nLL(betas, XN1f(itrain,:), XN1s(itrain,:), Xk(itrain,:), t(itrain), y(itrain)) + ...
                           lambda(i)*abs(betas)'*beta_scale;
            btrain = fmincon(fun, b0, A, b, [], [], [], [], [], opts);
            betas_est_lambda{i}(:,k) = btrain;
            nLL_est_lambda{i}(k) = calc_nLL(btrain, XN1f(itest,:), XN1s(itest,:), Xk(itest,:), t(itest), y(itest));
        end
    end
    toc
    S.betas_est_lambda = betas_est_lambda;
    S.nLL_est_lambda = nLL_est_lambda;
    
    % assemble a table including all the information for estimating lambda
    T_lambda.lambda = lambda(:);
    T_lambda.LL_per_trode = cellfun(@(x) -1 * mean(x)/sum(i_est_lambda), nLL_est_lambda);
    for i = 1:n_mdl
        [~,idx] = min(abs(nLL_est_lambda{i}-median(nLL_est_lambda{i})));
        T_lambda.betas(i,:) = betas_est_lambda{i}(:, idx)';
    end
    T_lambda = struct2table(T_lambda);
    S.T_lambda = T_lambda;
    
    % estimate the parameter coefficients
    [~,i_opt_lambda] = max(T_lambda.LL_per_trode);
    S.i_opt_lambda = i_opt_lambda;
    XN1f = X.N1f(~i_est_lambda,:);
    XN1s = X.N1s(~i_est_lambda,:);
    Xk = X.k(~i_est_lambda,:); 
    t = T_trode.days_since_init(~i_est_lambda);
    y = T_trode.(P_in.metric)(~i_est_lambda);
    fun = @(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y) + ...
                   lambda(i_opt_lambda)*abs(betas)'*beta_scale;
    betas = fmincon(fun, b0, A, b, [], [], [], [], [], opts);
    
    % collect variables into the output structure S
    S.P_in = P_in;
    S.T_trode = T_trode;
    S.T_regressor = T_regressor;
    S.T_dsgn = T_dsgn;
    S.T_mdl = T_mdl;
    S.betas = betas;
    
    S.T_betas = cell2table(num2cell(S.betas'));
    S.T_betas.Properties.VariableNames = ['kf', 'ks', P_in.model_parameters];
end
%% GENERATE_RESPONSE_VAR
%
%     
function [ygen, bgen] = generate_response_var(X, t)

    n_XN1f = size(X.N1f,2)-1;
    n_XN1s = size(X.N1s,2)-1;
    n_Xk = size(X.k,2); 

    B.kf = -rand;
    B.ks = -rand/100;
    B.N1f0 = rand;
    B.N1f = rand(n_XN1f,1);
    B.N1f(randperm(n_XN1f, floor(n_XN1f/2))) = 0;
    B.N1s0 = rand;
    B.N1s = rand(n_XN1s,1);
    B.N1s(randperm(n_XN1s, floor(n_XN1s/2))) = 0;
    B.k = -rand(n_Xk,1)/10;
    B.k(randperm(n_Xk, floor(n_Xk/2))) = 0;
    
    bgen = cell2mat(struct2cell(B));
    lambda = calc_resp_var(bgen, X.N1f, X.N1s, X.k, t);
    ygen = poissrnd(lambda);
    ygen = lambda;
    % add some noise
%     n = numel(ygen);                    
%     ygen(randperm(n,floor(n/4))) = poissrnd(mean(lambda), floor(n/4),1);
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