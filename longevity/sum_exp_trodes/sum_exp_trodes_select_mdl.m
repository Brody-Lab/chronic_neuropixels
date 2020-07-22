% SUM_EXP_TRODES_SELECT_MDL Best-subset variable selection using the
% sum-of-exponentials model fit to the unit count of each electrode in each
% recording.
%
%   N = Pois(lambda)
%
%   lambda = N1 * (alpha*exp(-(t-1)/tau_alpha) + (1-alpha)*exp(-(t-1)/tau))
%
%   tau = b0 + b_AP*AP + b_DV*DV + b_ML*ML + b_SVP&SVP + b_SPA*SPA;
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
%   iterations
%       The number of times the data are split in half, one part used for
%       model selection and another for parameter estimation.
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
%   shuffle
%       A logical scalar specifying whether to shuffle the rows for each
%       variable independently.
%       
function S = sum_exp_trodes_select_mdl(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'iterations', 10, ...
            @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'KFold', 5, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.default_model_parameters, ...
        @(x) all(ismember(x, P.sum_exp_trodes.possible_model_parameters)))
    addParameter(parseobj, 'noise', 'gaussian', @(x) any(strcmpi(x, {'poisson', 'gaussian'})))
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
    T_mdl = make_T_mdl('model_parameters', P_in.model_parameters);
    T_dsgn = make_design_matrix(T_trode, 'model_parameters', P_in.model_parameters);
    switch P_in.noise
        case 'gaussian'
            crossval_mdl = @crossval_gaussian;
            P_in.err_txt = 'MSE';
        case 'poisson'
            crossval_mdl = @crossval_poisson;
            P_in.err_txt = 'nLL';
        otherwise
            error('unrecognized noise distribiution');
    end
    % set up the (par)for loop
    % split the data and fit half of the data to the model variants
    Err = nan(size(T_mdl,1), P_in.iterations);
    n_trodes = size(T_trode,1);
    trodes_used = false(n_trodes, P_in.iterations);
    b = arrayfun(@(x) nan(numel(P_in.model_parameters)+3, P_in.iterations), (1:size(T_mdl,1))', 'uni', 0);
    for i = 1:P_in.iterations
        i_trodes=ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
        trodes_used(:,i)=i_trodes;
        y = T_trode.(P_in.metric)(i_trodes,:);
        t = T_trode.days_since_init(i_trodes,:);
        cvp = cvpartition(sum(i_trodes), 'KFold', P_in.KFold); % just to save a bit of time
        parfor j = 1:size(T_mdl,1)     
            [XN1, Xk] = partition_T_dsgn(T_dsgn, T_mdl, j, 'i_trodes', i_trodes);
            % fit the parameters using all the data in the subset
            betas = fit_mdl_sum_exp_trodes(XN1, Xk, t, y, 'noise', P_in.noise);
            b{j}(1:3,i) = betas(1:3);
            betas_regressors = nan(numel(T_mdl{j,:}),1);
            betas_regressors(T_mdl{j,:}) = betas(4:end);
            b{j}(4:end,i) = betas_regressors;
            % get out of sample LL
            err = crossval(crossval_mdl, XN1, Xk, t, y, 'Partition', cvp);
            Err(j,i) = mean(err);
            fprintf('\nIteration %i - %s model %i - %s: %0.3f', ...
                    i, P_in.noise, j, P_in.err_txt, Err(j,i))
        end
    end
    % collect variables into the table T_MDL
    % in each iteration, get the rank of that model from the lowest to
    % highest LL
    Err=real(Err);
    err_norm = (Err - min(Err))./(max(Err)-min(Err));
    T_res = T_mdl;
    T_res.err_norm = mean(err_norm,2);
    T_res.err_norm(isnan(T_res.err_norm)) = 0;
    T_res.b = b;
    T_res.MSE = Err;    
    [~, mdl_sort_idx] = sort(Err); 
    for i = 1:P_in.iterations
        T_res.rank(mdl_sort_idx(:,i),i) = (1:size(err_norm,1))';
    end
    T_res.frac_best = sum(T_res.rank==1,2)/P_in.iterations;
    T_res.avg_rank = mean(T_res.rank,2);
    T_res.n_regressors = sum(T_res{:,P_in.model_parameters},2);
    
    % sort rows
    sort_col = find(strcmp(T_res.Properties.VariableNames,'err_norm'));
    [T_res, I] = sortrows(T_res, sort_col); 
    T_mdl = T_mdl(I, :);
    
    % collect variables into the output structure S
    S.T_trode = T_trode;
    S.T_regressor = T_regressor;
    S.trodes_used = trodes_used;
    S.T_res = T_res;
    S.T_mdl = T_mdl;
    S.P_in = P_in;
    S.T_dsgn = T_dsgn;
end
%% CROSSVAL_gaussian
function mse = crossval_gaussian(XN1train, Xktrain, ttrain,  ytrain, ...
                               XN1test, Xktest,  ttest, ytest)
    b = fit_mdl_sum_exp_trodes(XN1train, Xktrain, ttrain, ytrain);
    yhat = predict_y(b, XN1test, Xktest, ttest);
    mse = mean((ytest-yhat).^2);
end

%% CROSSVAL_poisson
function nLL = crossval_poisson(XN1train, Xktrain, ttrain,  ytrain, ...
                               XN1test, Xktest,  ttest, ytest)
    b = fit_mdl_sum_exp_trodes(XN1train, Xktrain, ttrain, ytrain);
    yhat = predict_y(b, XN1test, Xktest, ttest);
    nLL = sum(yhat) - sum(ytest.*log(yhat)); % ignoring the term independent of lambda
end