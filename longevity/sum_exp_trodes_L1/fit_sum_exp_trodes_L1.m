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
%   shuffle
%       A logical scalar specifying whether to shuffle the rows for each
%       variable independently.
%       
function S = fit_sum_exp_trodes_L1(Cells, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 10, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.model_parameters_L1, ...
        @(x) all(ismember(x, P.sum_exp_trodes.model_parameters_L1)))
    addParameter(parseobj, 'noise', 'poisson', @(x) any(strcmpi(x, {'poisson', 'gaussian'})))
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
    % fit the regression without regularization to estimate the relative
    % scales of the penalization parameters for each each model parameter
    B0 = fit_SoE_mdl(X, T_trode.days_since_init, T_trode.(P_in.metric), 'noise', P_in.noise);
    penal_N1f = sum(abs(B0.N1f));
    penal_N1s = sum(abs(B0.N1s));
    penal_N1k = sum(abs(B0.N1k));
    
    [B,FitInfo] = fit_SoE_mdl_L1(X, T_trode.days_since_init, T_trode.(P_in.metric), ...
                       'noise', P_in.noise);
    
    % collect variables into the output structure S
    S.P_in = P_in;
    S.T_trode = T_trode;
    S.T_regressor = T_regressor;
    S.T_dsgn = T_dsgn;
    S.T_mdl = T_mdl;
    S.B = B;
    S.FitInfo = FitInfo;
end