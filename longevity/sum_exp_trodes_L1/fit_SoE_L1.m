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
%   cvp
%       CV partition
%
%   generate_data
%       Make fake response data using the observed predictors
%
%   infer_betas
%       A scalar toggle for whether to infer the model coefficients (as
%       opposed to just returning the optimal lambda, in-sample LL, and
%       out-of-sample LL)
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
%       penalization paer_ramter
%
%   shuffle
%       A logical scalar specifying whether to shuffle the rows for each
%       variable independently.
%       
function S = fit_SoE_L1(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'cvp', [])
    addParameter(parseobj, 'infer_betas', true, @(x) isscalar(x) && (x == 0 || x== 1))
    addParameter(parseobj, 'KFold', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'generate_data', false, @(x) x==0 || x==1)
    addParameter(parseobj, 'i_est_lambda', [])
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.params_N1f_N1s, ...
        @(x) all(ismember(x, P.sum_exp_trodes.possible_model_parameters)))
    addParameter(parseobj, 'noise', 'poisson', @(x) any(strcmp(x, {'poisson', 'gaussian'})))
    addParameter(parseobj, 'num_lambda', 12, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}))
    addParameter(parseobj, 'lambda_ratio', 1e-2, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}))
    addParameter(parseobj, 'shuffle', false, @(x) x==0 || x==1)
    addParameter(parseobj, 'T_trode', [], @(x) istable(x))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    S.P_in = P_in;
    param_types = arrayfun(@get_parameter_type, P_in.model_parameters, 'uni', 0);
    if ismember("N1", param_types) && any(ismember(["N1s", "N1f"], param_types))
        error('A model cannot have both the N1 term and either the N1s or N1f term.')
    end
    if isempty(P_in.T_trode)
        [T_trode, T_regressor] = ...
                  make_T_trode(Cells, 'model_parameters', P_in.model_parameters, ...
                                      'normalize_regressors', true, ...
                                      'unit_distance', P.unit_distance, ...
                                      'x0', P.x0);
        S.T_regressor = T_regressor;
    else
        T_trode = P_in.T_trode;
    end
    
    if P_in.shuffle
        n = size(T_trode,1);
        for i = 1:numel(T_regressor.name)
            T_trode.(T_regressor.name{i}) = T_trode.(T_regressor.name{i})(randperm(n));
        end
    end
    
    SX = make_separate_design_matrices(T_trode, P_in.model_parameters);
    S.SX = SX;
    
    if P_in.generate_data
        [ygen, bgen] = generate_response_var(SX,  T_trode.days_since_init);
        T_trode.(P_in.metric) = ygen;
        S.bgen = bgen;
    end
    S.T_trode = T_trode;
    
    
    switch P_in.noise
        case 'gaussian'
            if isfield(SX, 'N1')
                calc_nLL = @calc_nLL_gaussian_N1_a;
            else
                calc_nLL = @calc_nLL_gaussian_N1f_N1s;
            end
        case 'poisson'
            if isfield(SX, 'N1')
                calc_nLL = @calc_nLL_poisson_N1_a;
            else
                calc_nLL = @calc_nLL_poisson_N1f_N1s;
            end
    end
    
    b0 = get_initial_estimates(SX);
    [A,b] = get_linear_inequality_constraints(SX);
          
    % randomly select half of the data for estimating lambda
    n_trodes = size(SX.k,1);
    if isempty(P_in.i_est_lambda)
        i_est_lambda = ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
    else
        i_est_lambda = P_in.i_est_lambda;
    end
    S.i_est_lambda = i_est_lambda;
    t = T_trode.days_since_init(i_est_lambda);
    y = T_trode.(P_in.metric)(i_est_lambda);
    Xk = SX.k{i_est_lambda,:}; 

    opts=optimoptions(optimoptions('fmincon'), ...
                        'Display', 'notify', ...
                        'MaxFunctionEvaluations', 1e5);
    
    %fit without penalization
    if isfield(SX, 'N1s')
        XN1f = SX.N1f{i_est_lambda,:};
        XN1s = SX.N1s{i_est_lambda,:};        
        XN1f = padarray(XN1f, [0,1], 1, 'pre');
        XN1s = padarray(XN1s, [0,1], 1, 'pre');
        
        [beta_no_reg, nLL_no_reg] = ...
            fmincon(@(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y), b0, A, b, ...
                                      [], [], [], [], [], ...
                                      opts);
    else
        XN1 = SX.N1{i_est_lambda,:};
        XN1 = padarray(XN1, [0,1], 1, 'pre');
        
        [beta_no_reg, nLL_no_reg] = ...
            fmincon(@(betas) calc_nLL(betas, XN1, Xk, t, y), b0, A, b, ...
                                      [], [], [], [], [], ...
                                      opts);
    end
    assert(isreal(nLL_no_reg))
    
    beta_weight = get_coefficient_weights(SX, beta_no_reg);
    S.beta_weight = beta_weight;
            
    loghi = -4; % this is empirical
    loglo = loghi+log(P_in.lambda_ratio);
    lambda = exp(linspace(loghi,loglo,P_in.num_lambda-1))*nLL_no_reg/(abs(beta_no_reg')*beta_weight);
    lambda = [lambda,0];
     
    % fit with penalization
    n_mdl = size(lambda,2);
    if isempty(P_in.cvp)
        cvp = cvpartition(sum(i_est_lambda), 'KFold', P_in.KFold);
    else
        cvp = P_in.cvp;
    end
    betas_est_lambda = cell(n_mdl,1);
    LL_in_sample = cell(n_mdl, 1);
    LL_per_trode = cell(n_mdl, 1);
    if isfield(SX, 'N1s')
        bundled_inputs = {XN1f, XN1s, Xk, t, y};
    else
        bundled_inputs = {XN1, Xk, t, y};
    end
    parfor i = 1:n_mdl
        betas_est_lambda{i} = nan(numel(b0), P_in.KFold);
        % To faciliate debugging during cross-validation, CROSSVAL is not used
        for k = 1:P_in.KFold
            itrain = training(cvp,k);
            itest = test(cvp,k);
            training_inputs = cellfun(@(x) x(itrain,:), bundled_inputs, 'uni', 0);
            testing_inputs = cellfun(@(x) x(itest,:), bundled_inputs, 'uni', 0);
            fun = @(betas) calc_nLL(betas, training_inputs{:}) + lambda(i)*abs(betas)'*beta_weight;
            btrain = fmincon(fun, b0, A, b, [], [], [], [], [], opts);
            betas_est_lambda{i}(:,k) = btrain;
            % out of sample loglikelihood
            LL_per_trode{i}(k) = calc_full_LL(P_in.noise, btrain, testing_inputs{:})/sum(itest);
            % in-sample loglikliehood, including the constant term, for calculating BIC
            LL_in_sample{i}(k,1) = calc_full_LL(P_in.noise, btrain, training_inputs{:});
        end
    end
    % assemble a table including all the information for estimating lambda
    T_lambda.lambda = lambda(:);
    T_lambda.LL_per_trode = cellfun(@mean, LL_per_trode);
        % the averaging is done across cv-folds
    assert(isreal(T_lambda.LL_per_trode))
    T_lambda.LL_in_sample = cellfun(@mean, LL_in_sample);
    for i = 1:n_mdl
        [~,idx] = min(abs(LL_per_trode{i}-median(LL_per_trode{i})));
        T_lambda.betas(i,:) = betas_est_lambda{i}(:, idx)';
    end
    T_lambda = struct2table(T_lambda);
    S.T_lambda = T_lambda;
    [~,i_opt_lambda] = max(T_lambda.LL_per_trode);
    S.i_opt_lambda = i_opt_lambda;
    % calculate BIC
    ntrain = arrayfun(@(x) sum(training(cvp,x)), (1:P_in.KFold)');
    S.BIC = mean(numel(b0)*log(ntrain) - 2*LL_in_sample{i_opt_lambda}); % averaging across folds
    %% estimate the parameter coefficients
    if ~P_in.infer_betas
        return
    end
    
    Xk = SX.k{~i_est_lambda,:}; 
    t = T_trode.days_since_init(~i_est_lambda);
    y = T_trode.(P_in.metric)(~i_est_lambda);

    if isfield(SX, 'N1s')
        XN1f = SX.N1f{~i_est_lambda,:};
        XN1s = SX.N1s{~i_est_lambda,:};
        XN1f = padarray(XN1f, [0,1], 1, 'pre');
        XN1s = padarray(XN1s, [0,1], 1, 'pre');
        fun = @(betas) calc_nLL(betas, XN1f, XN1s, Xk, t, y) + ...
                       lambda(i_opt_lambda)*abs(betas)'*beta_weight;
    else
        XN1 = SX.N1{~i_est_lambda,:};
        XN1 = padarray(XN1, [0,1], 1, 'pre');
        fun = @(betas) calc_nLL(betas, XN1, Xk, t, y) + ...
                       lambda(i_opt_lambda)*abs(betas)'*beta_weight;
    end           
    betas = fmincon(fun, b0, A, b, [], [], [], [], [], opts);
    S.betas = betas;
    
    % make a table for the betas so that they can be more easily readout
    S.T_betas = cell2table(num2cell(S.betas'));
    if isfield(SX, 'N1s')
        S.T_betas.Properties.VariableNames = ['kf', 'ks', 'N1f0', 'N1s0', P_in.model_parameters];
    else
        S.T_betas.Properties.VariableNames = ['kf', 'ks', 'a', 'N1', P_in.model_parameters];
    end
end
%% CALC_FULL_LL
%
%   Calculate the summed loglikelhood including the constant terms
%
%=INPUT, Positional
%
%   1) noise
%       A char array indicating whether the noise distribution is Poisson
%       or Gaussian (either "Poisson" or "Gaussian", case-insensitive)
%
%   2-7) If there are seven inputs total, then they are assumed to be:
%       2) b - coefficients
%       3) XN1f - design matrix for regressors in the N1f term
%       4) XN1s - design matrix for regressors in the N1s term
%       5) Xk - design matrix for regressors in the k term
%       6) t - a vector of days elapsed
%       7) y - the response variable, a vector
%
%   2-6) If there are six inputs in total, then they are assumed to be:
%       2) b - coefficients
%       3) XN1 - design matrix for regressors in the N1 term
%       4) Xk - design matrix for regressors in the k term
%       5) t - a vector of days elapsed
%       6) y - the response variable, a vector
%
%=OUTPUT
%
%   LL
%       A scalar indicating the loglikelihood
%
%=EXAMPLE CALL
%
%   >>LL = calc_full_LL('poisson', b, XN1f, XN1s, Xk, t, y)
%   >>LL = calc_full_LL('gaussian', b, XN1, Xk, t, y)
function LL = calc_full_LL(noise, varargin)
    assert(strcmpi(noise, 'poisson') || strcmpi(noise, 'gaussian'))
    if numel(varargin) == 6
        b = varargin{1};
        XN1f = varargin{2};
        XN1s = varargin{3};
        Xk = varargin{4};
        t = varargin{5};
        y = varargin{6};
        n = numel(t);
        if strcmp(noise, 'poisson')
            LL = -calc_nLL_poisson_N1f_N1s(b, XN1f, XN1s, Xk, t, y) - mean(log(factorial(y)));
        else
            yhat = calc_resp_var_N1f_N1s(b, XN1f, XN1s, Xk, t);
            sigma2hat = 1/(-2) * sum((y - yhat).^2);
            sq_err = calc_nLL_gaussian_N1f_N1s(b, XN1f, XN1s, Xk, t, y);
            LL = -0.5*numel(t)*(log(sigma2hat) + log(2*pi)) - ...
                 1/2/sigma2hat*sq_err;
        end
    elseif numel(varargin) == 5
        b = varargin{1};
        XN1 = varargin{2};
        Xk = varargin{3};
        t = varargin{4};
        y = varargin{5};
        n = numel(t);
        if strcmp(noise, 'poisson')
            LL = -calc_nLL_poisson_N1_a(b, XN1, Xk, t, y) - mean(log(factorial(y)));
        else
            yhat = calc_resp_var_N1_a(b, XN1, Xk, t);
            sigma2hat = 1/(-2) * sum((y - yhat).^2);
            sq_err = calc_nLL_gaussian_N1_a(b, XN1, Xk, t, y);
            LL = -0.5*numel(t)*(log(sigma2hat) + log(2*pi)) - ...
                 1/2/sigma2hat*sq_err;
        end
    else
        error('Cannot parse inputs')
    end
end
%% GENERATE_RESPONSE_VAR
%   Simulate the response variable using random model coefficients and the
%   observed predictors
%     
%=INPUT
%
%   SX
%       A structure of design matrices for each term in the model
%
%   t
%       Number of days from implantation minus one
%
%=OUTPUT
%
%   ygen
%       generated response variable
%
%   bgen
%       generated model coefficients
%
function [ygen, bgen] = generate_response_var(SX, t, varargin)
    parseobj = inputParser;
    addParameter(parseobj, 'no_noise', false, @(x) x==0 || x==1)
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
 
    if all(isfield(SX, {'N1f', 'N1s', 'k'})) && ...
       ~isfield(SX, 'N1')
   
        B.kf = -rand;
        B.ks = -rand/100;
        B.N1f0 = rand;
        B.N1s0 = rand;
        
        n_XN1f = size(SX.N1f,2)-1;
        n_XN1s = size(SX.N1s,2)-1;
        
        B.N1f = rand(n_XN1f,1);
        B.N1f(randperm(n_XN1f, floor(n_XN1f/2))) = 0;
        B.N1s = rand(n_XN1s,1);
        B.N1s(randperm(n_XN1s, floor(n_XN1s/2))) = 0;
        
        n_Xk = size(X.k,2); 
        B.k = -rand(n_Xk,1)/10;
        B.k(randperm(n_Xk, floor(n_Xk/2))) = 0;

        bgen = cell2mat(struct2cell(B));
        lambda = calc_resp_var_N1f_N1s(bgen, SX.N1f, SX.N1s, SX.k, t);
        
    elseif all(isfield(SX, {'N1', 'k'})) && ...
           all(~isfield(SX, {'N1f', 'N1s'}))

        B.kf = -rand;
        B.ks = -rand/100;
        B.a = rand;
        B.N10 = rand;
        
        n_XN1 = size(SX.N1,2)-1;
        B.N1 = rand(n_XN1,1);
        B.N1(randperm(n_XN1, floor(n_XN1/2))) = 0;
        
        n_Xk = size(X.k,2); 
        B.k = -rand(n_Xk,1)/10;
        B.k(randperm(n_Xk, floor(n_Xk/2))) = 0;

        bgen = cell2mat(struct2cell(B));
        lambda = calc_resp_var_N1_a(bgen, SX.N1f, SX.N1s, SX.k, t);
    else
        error('Unknown structure of design matrices')
    end
    
    if P_in.no_noise
        ygen = lambda;
    else
        ygen = poissrnd(lambda);
    end
end
%% GET_COEFFICIENT_WEIGHTS
%   return a vector of weights for each coefficient used to calculate the
%   penalization term
%
%=INPUT
%
%   SX
%       A structure of design matrices for each term in the model
%
%   beta_no_reg
%       A vector of coefficient estimates made without regularization
%
%=OUTPUT
%
%   wts
%       A vector of weights with the same number of elements as the number
%       of coefficients

function wts = get_coefficient_weights(SX, beta_no_reg)

    if all(isfield(SX, {'N1f', 'N1s', 'k'})) && ...
       ~isfield(SX, 'N1')
   
        n_N1f = size(SX.N1f,2);
        n_N1s = size(SX.N1s,2);
        n_k   = size(SX.k,2);
        B.kf   = beta_no_reg(1);
        B.ks   = beta_no_reg(2);
        B.N1f0 = beta_no_reg(3);
        B.N1s0 = beta_no_reg(4);
        B.N1f  = beta_no_reg(5:n_N1f+4);
        B.N1s  = beta_no_reg(n_N1f+5:n_N1f+n_N1s+4);
        B.k    = beta_no_reg(n_N1f+n_N1s+5:end);
        assert(numel(cell2mat(struct2cell(B))) == numel(beta_no_reg));
        s = mean(abs([B.N1f; B.N1s]))/mean(abs(B.k));
        wts = [0; 0; 0; 0;
               ones(n_N1f,1); ...
               ones(n_N1s,1); ...
               s*ones(n_k,1)];
        
    elseif all(isfield(SX, {'N1', 'k'})) && ...
           all(~isfield(SX, {'N1f', 'N1s'}))
       
        n_N1 = size(SX.N1,2);
        n_k  = size(SX.k,2);
        B.kf  = beta_no_reg(1);
        B.ks  = beta_no_reg(2);
        B.a   = beta_no_reg(3);
        B.N10 = beta_no_reg(4);
        B.N1  = beta_no_reg(5:n_N1+4);
        B.k   = beta_no_reg(n_N1+5:end);
        assert(numel(cell2mat(struct2cell(B))) == numel(beta_no_reg));
        s = mean(abs(B.N1))/mean(abs(B.k));
        wts = [0; 0; 0; 0;
               ones(n_N1,1); ...
               s*ones(n_k,1)];
    else
        error('Unknown structure of design matrices')
    end
end
%% GET_INITIAL_ESTIMATES
%   Return the initial estimates of the coefficients
%
%=INPUT
%
%   SX
%       A structure of design matrices for each term in the model
%
%=OUTPUT
%
%   b0
%       A vector of initial estimates of the parameter values
function b0 = get_initial_estimates(SX)
    if all(isfield(SX, {'N1f', 'N1s', 'k'})) && ...
       ~isfield(SX, 'N1')
   
        n_N1f = size(SX.N1f,2);
        n_N1s = size(SX.N1s,2);
        n_k   = size(SX.k,2);
        b0 = [-1; -0.01; 1; 1; ones(n_N1f,1); ones(n_N1s,1); -0.01*ones(n_k,1)]; 
        % b0 = [kf; ks; N1f0; N1s0; ...
    elseif all(isfield(SX, {'N1', 'k'})) && ...
           all(~isfield(SX, {'N1f', 'N1s'}))

        n_N1 = size(SX.N1,2);
        n_k   = size(SX.k,2);
        b0 = [-1; -0.01; 0.5; 1; ones(n_N1,1); -0.01*ones(n_k,1)];
        % b0 = [kf; ks; a; N1; ...
    else
        error('Unknown structure of design matrices')
    end
end
%% GET_LINEAR_INEQUALITY_CONSTRAINTS
%   Return the linear inequality constraints expressed as an M-by-N matrix
%   A and an M-element vector b.
%
%=INPUT
%
%   SX
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
function [A,b] = get_linear_inequality_constraints(SX)

    if all(isfield(SX, {'N1f', 'N1s', 'k'})) && ...
       ~isfield(SX, 'N1')
   
        n_N1f = size(SX.N1f,2);
        n_N1s = size(SX.N1s,2);
        n_k   = size(SX.k,2);

        % make k_fast <= k_slow
        A_kfs = [1,-1]; % k_fast - k_slow <= 0

        % make sure the N1s and N1f terms are positive for all regressors (who
        % are in the range of [0,1]). Negative because the constraint's format
        % is A*x <= b, and we want the linear combination to be >= 0
        if n_N1f > 1
            A_N1f = -double(dec2bin(0:2^(n_N1f)-1)=='1');
            A_N1f = [-ones(size(A_N1f,1),1), A_N1f]; % the constant term
        else
            A_N1f = -1; % the constant term
        end
        if n_N1s > 1
            A_N1s = -double(dec2bin(0:2^(n_N1s)-1)=='1');
            A_N1s = [-ones(size(A_N1s,1),1), A_N1s]; % the constant term
        else
            A_N1s = -1;
        end
        A = blkdiag(A_kfs, A_N1f, A_N1s);
        A = padarray(A, [0,n_k], 0, 'post');
        b = zeros(size(A,1),1);
        
        % rearrange the columns so that they correspond to kf, ks, N1f0,
        % N1s0, ....
        iN1s0 = size(A_kfs,2)+size(A_N1f,2)+1;
        A = [A(:,1:3), A(:,iN1s0), A(:, 4:(iN1s0-1)), A(:,iN1s0+1:end)];
        
    elseif all(isfield(SX, {'N1', 'k'})) && ...
           all(~isfield(SX, {'N1f', 'N1s'}))

        n_N1 = size(SX.N1,2);
        n_k   = size(SX.k,2);

        % make k_fast <= k_slow
        A_kfs = [1,-1]; % k_fast - k_slow <= 0
        
        if n_N1 > 1
            A_N1 = -double(dec2bin(0:2^(n_N1)-1)=='1');
            A_N1 = [-ones(size(A_N1,1),1), A_N1]; % the constant term
        else
            A_N1 = -1; % the constant term
        end
        A = blkdiag([A_kfs, 0], A_N1);
        A = padarray(A, [0,n_k], 0, 'post');
        b = zeros(size(A,1),1);
    else
        error('Unknown structure of design matrices')
    end
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