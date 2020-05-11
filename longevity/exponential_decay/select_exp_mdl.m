% SELECT_EXP_MDL Variable selection in the exponential model using best
% subset selection

%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILE and POSTPROCESS_CELLS
%
%=OUTPUT
%
%   S
%       A structure with the following fields
%       - T_mdl, a table of models
%       - idx_mdl_selection, the index of electrode recordings used for
%       model selection. The electrode-recordings not used for model
%       selection are used for parameter estimation
%       - b_hat, the estimated coefficients from the optimal model, done
%       without bootstrapping
%       - b_boot, a matrix of coefficients, in which each row is a parameter
%       estimate and each column is an iteration of bootstrap sampling from
%       the data. Parameter estimation was performed for the only the model
%       selected in the model selection stage. 
%       - b_names, the array of regressor names for the selected model
%       - b_ci, the Studentized bootstrap confidence intervals, 
%       
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   iterations
%       The number of times the model is split
%
%   KFold
%       Number of cross-validation folds
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   n_boot
%       A positive integer specifying the number of iterations of bootstrap
%       draws.
%
%   unit_distance
%       A scalar nonnegative, integer indicating the maximum distance (in
%       electrode index) between an electrode and the site where the peak
%       amplitude of a unit was observed for that unit to be counted for
%       that electrode. A value 0 indicates that a unit is counted for an
%       electrode only if that unit's peak response occurred on that
%       electrode.
%
%   x0
%       The first day after surgery to be examined.
%
function S = select_exp_mdl(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 5, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'iterations', 10, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    T_trode = make_T_trode(Cells, 'x0', P.x0, ...
                                  'unit_distance', P.unit_distance, ...
                                  'mean_subtract_factors', true, ...
                                  'subtract_x0', true);
    factor_range = get_range_of_trode_factor(T_trode);
    T_trode{:,P.ED_trode_factors}=T_trode{:,P.ED_trode_factors}./...
                                 max(T_trode{:,P.ED_trode_factors});
    y = T_trode.(P_in.metric);                                
    T_mdl = select_ED_model_subsets;                                    
    % set up the (par)for loop
    % split the data and fit half of the data to the model variants
    mse = nan(size(T_mdl,1), P_in.iterations);
    n_trodes = size(T_trode,1);
    trodes_used = false(n_trodes, P_in.iterations);
    for i = 1:P_in.iterations
        trodes_used(:,i)=ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
        yi = y(trodes_used(:,i),:);
        t = T_trode.days_elapsed(trodes_used(:,i));
        cvp = cvpartition(sum(trodes_used(:,i)), 'KFold', P_in.KFold); % just to save a bit of time
        parfor j = 1:size(T_mdl,1)
            [X0, Xt] = make_X0_Xt(T_trode, T_mdl{j,:}, trodes_used(:,i));
            mse(j,i) = crossval('mse', X0, Xt, t, yi, ...
                                'PredFun', @fit_and_pred_ed_fmincon,  ...
                                'Partition', cvp);
            fprintf('\nIteration %i - model %i - mse: %0.3f', i, j, mse(j,i))
        end
    end
    % collect variables into the table T_MDL
    % in each iteration, get the rank of that model from the lowest to
    % highest mse
    T_mdl.mse = mse;    
    [~, mdl_sort_idx] = sort(mse);
    for i = 1:P_in.iterations
        T_mdl.rank(mdl_sort_idx(:,i),i) = (1:size(mse,1))';
    end
    T_mdl.frac_best = sum(T_mdl.rank==1,2)/P_in.iterations;
    T_mdl.avg_rank = mean(T_mdl.rank,2);
    T_mdl.mse_med = median(mse,2);
    T_mdl.n_regressors = sum(T_mdl{:,P.ED_trode_regressors_all},2);
    
    % sort rows
    sort_col = find(strcmp(T_mdl.Properties.VariableNames,P.ED_trode_selection_criterion));
    T_mdl = sortrows(T_mdl, sort_col); 
    
    % collect variables into the output structure S
    S.T_trode = T_trode;
    S.trodes_used = trodes_used;
    S.T_mdl = T_mdl;
    S.P_in = P_in;
    S.factor_range = factor_range;
end
%% FIT_AND_PRED_ED_FMINCON
% fit an exponential function using FMINUNC and predict values
%
%=INPUT
%
%   Xtrain
%       The columns of the design matrix with values for the regressors,
%       for training.
%
%   tDXtrain
%       The columns of the design that are the product of values of the
%       regressors multiplied to the numbers of days elapsed, for training.
%
%   ytrain
%       The observed unit count, for training.
%
%   Xtest
%       The columns of the design matrix with values for the regressors,
%       for testing.
%
%   tDXtest
%       The columns of the design that are the product of values of the
%       regressors multiplied to the numbers of days elapsed, for testing.
%
%=OUTPUT
%   
%   yhat
%       values predicted by XTEST and XTTEST
function yhat = fit_and_pred_ed_fmincon(X0train, Xttrain, ttrain, ytrain, X0test, Xttest, ttest)
    opts=optimset('fmincon');
    opts.TolFun = 1e-10;
    opts.Display = 'off';
    % ones for the regressors associated with initial count (y0) and 0's for
    % the regressors associated with the decay rate (k)
    n_X = size(X0train,2);
    n_Xt = size(Xttrain,2);
    b0 = [ones(n_X,1); ones(n_Xt,1)]; 
    % A*b >= 0
    % The weighted sum of the decay terms has to be positive.
    % All values of X are normalized to be in [0,1]
    if n_Xt > 1
        A_Xt = double(dec2bin(0:2^(n_Xt-1)-1)=='1');
        A_Xt = [ones(size(A_Xt,1),1), A_Xt];
    else
        A_Xt = 1;
    end
    ncol_A = size(A_Xt,1);
    A = [zeros(ncol_A,n_X), -A_Xt];
    b = fmincon(@(b) ls_exp(b, X0train, Xttrain, ttrain, ytrain), b0, A, zeros(ncol_A,1), ...
                            [], [], [], [], [], ...
                            opts);
    yhat = eval_exp(b, X0test, Xttest, ttest);
%     debug_exp_fit(b,X0train,Xttrain,ttrain, ytrain)
%     n = size(X0train,2);
%     median(Xttest*b(n+1:end))
end
%% Caclulate the sum of least squares
function Q = ls_exp(b,X0,Xt,t, y)
    yhat = eval_exp(b,X0,Xt, t);
    Q = sum((y - yhat).^2);
end
%% Evaluate a exponential function
function yhat = eval_exp(b,X0,Xt,t)
    n = size(X0,2);
    y0 = X0*b(1:n);
    tDtau = t ./ (Xt*b(n+1:end));
    yhat = y0.*exp(-tDtau);   
end