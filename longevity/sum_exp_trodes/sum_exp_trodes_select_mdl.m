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
%       - T_mdl, a table of models
%       - T_trode, a table of each electrode and each recoding
%       - factor_range, the range of each experimental factor in T_trode,
%       before being normalized to between [0,1] in T_trode. This is
%       required to knowing the coefficient weight of an experimental
%       factor per the unit of that factor.
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
function S = sum_exp_trodes_select_mdl(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 5, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    addParameter(parseobj, 'iterations', 10, ...
            @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'shuffle', false, @(x) x==0 || x==1)
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    [T_trode, exp_factors] = ...
              make_T_trode(Cells, 'x0', P.x0, ...
                                  'unit_distance', P.unit_distance, ...
                                  'normalize_factors', true);
    if P_in.shuffle
        n = size(T_trode,1);
        for i = 1:numel(exp_factors.name)
            T_trode.(exp_factors.name{i}) = T_trode.(exp_factors.name{i})(randperm(n));
        end
    end
    T_mdl = make_T_mdl; % that needs to be changed to include the constant terms
    T_dsgn = make_design_matrix(T_trode);
    
    % set up the (par)for loop
    % split the data and fit half of the data to the model variants
    MSE = nan(size(T_mdl,1), P_in.iterations);
    n_trodes = size(T_trode,1);
    trodes_used = false(n_trodes, P_in.iterations);
    b = arrayfun(@(x) nan(numel(P.sum_exp_trodes.regressors)+2, P_in.iterations), (1:size(T_mdl,1))', 'uni', 0);
    for i = 1:P_in.iterations
        i_trodes=ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
        trodes_used(:,i)=i_trodes;
        y = T_trode.(P_in.metric)(i_trodes,:);
        t = T_trode.days_since_init(i_trodes,:);
        cvp = cvpartition(sum(i_trodes), 'KFold', P_in.KFold); % just to save a bit of time
        parfor j = 1:size(T_mdl,1)     
            [XN1, Xk] = partition_T_dsgn(T_dsgn, T_mdl, j, 'i_trodes', i_trodes);
            % fit the parameters using all the data in the subset
            betas = fit_mdl_sum_exp_trodes(XN1, Xk, t, y);
            b{j}(1:2,i) = betas(1:2);
            betas_regressors = nan(numel(T_mdl{j,:}),1);
            betas_regressors(T_mdl{j,:}) = betas(3:end);
            b{j}(3:end,i) = betas_regressors;
            % get out of sample LL
            mse = crossval(@crossval_mdl, XN1, Xk, t, y, 'Partition', cvp);
            MSE(j,i) = mean(mse);
            fprintf('\nIteration %i - model %i - MSE: %0.3f', i, j, MSE(j,i))
        end
    end
    % collect variables into the table T_MDL
    % in each iteration, get the rank of that model from the lowest to
    % highest LL
    
    MSE_norm = (MSE - min(MSE))./(max(MSE)-min(MSE));
    T_res = T_mdl;
    T_res.MSE_norm = mean(MSE_norm,2);
    T_res.MSE_norm(isnan(T_res.MSE_norm)) = 0;
    T_res.b = b;
    T_res.MSE = MSE;    
    [~, mdl_sort_idx] = sort(MSE); 
    for i = 1:P_in.iterations
        T_res.rank(mdl_sort_idx(:,i),i) = (1:size(MSE,1))';
    end
    T_res.frac_best = sum(T_res.rank==1,2)/P_in.iterations;
    T_res.avg_rank = mean(T_res.rank,2);
    T_res.n_regressors = sum(T_res{:,P.sum_exp_trodes.regressors},2);
    
    % sort rows
    sort_col = find(strcmp(T_res.Properties.VariableNames,P.sum_exp_trodes.selection_criterion));
    [T_res, I] = sortrows(T_res, sort_col); 
    T_mdl = T_mdl(I, :);
    
    % collect variables into the output structure S
    S.T_trode = T_trode;
    S.exp_factors=exp_factors;
    S.trodes_used = trodes_used;
    S.T_res = T_res;
    S.T_mdl = T_mdl;
    S.P_in = P_in;
    S.T_dsgn = T_dsgn;
end
%% CROSSVAL_MDL
function mse = crossval_mdl(XN1train, Xktrain, ttrain,  ytrain, ...
                               XN1test, Xktest,  ttest, ytest)
    b = fit_mdl_sum_exp_trodes(XN1train, Xktrain, ttrain, ytrain);
    yhat = predict_y(b, XN1test, Xktest, ttest);
    mse = mean((ytest-yhat).^2);
end