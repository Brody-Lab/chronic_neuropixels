% est_coeff_of_ED_mdl estimate the coefficients of the regressors of the
% models specified in S.T_mdl
%
%=INPUT
%
%   S
%       A structure made using SELECT_EXP_MDL
%
%=OUTPUT
%
%   b
%       the coefficients of the selected model
%
%   i_mdl
%       A linear indexing array to specify which models for
%       fitting
%
%=OPTIONAL INPUT
function S = est_coeff_of_ED_mdl(S, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'n_boot',P.exp_decay_n_boots, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
    addParameter(parseobj, 'i_mdl',1, ...
        @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    y = S.T_trode.(S.P_in.metric);
    n_iter = size(S.trodes_used,2); % number of times that data was randomly split
    n_regressors = numel(P.ED_trode_regressors_all);
    n_mdls_tested = size(S.T_mdl,1);
    
    % add field
    for f = {'b', 'cil', 'ciu'}
        f = f{:};
        if ~isvar(S.T_mdl, f)
            S.T_mdl.(f) = nan(n_mdls_tested, n_regressors);
        end
    end
    for m = 1:numel(P_in.i_mdl)
        ind_m = P_in.i_mdl(m);
        idx_regressors = S.T_mdl{ind_m,1:n_regressors};
        b = nan(n_regressors, n_iter);
        cil = nan(n_regressors, n_iter);
        ciu = nan(n_regressors, n_iter);
        boot_coeff = cell(1,n_iter);
        for i = 1:n_iter
            yi = y(S.trodes_used(:,i),:);
            t = S.T_trode.days_elapsed(S.trodes_used(:,i));
            [X0, Xt] = make_X0_Xt(S.T_trode, idx_regressors, S.trodes_used(:,i));
            b(idx_regressors,i) = fit_ed(X0, Xt, t, yi);
            bootcoeff = nan(sum(idx_regressors), P_in.n_boot);
            n_trodes = sum(S.trodes_used(:,i));
            parfor j = 1:P_in.n_boot
                bootidx = datasample(1:n_trodes, n_trodes);
                bootcoeff(:,j) = fit_ed(X0(bootidx,:), Xt(bootidx,:),t(bootidx,:), yi(bootidx,:));
                fprintf('\n Iteration %i boot %i',i, j)
            end
            boot_coeff{i} = nan(n_regressors, P_in.n_boot);
            boot_coeff{i}(idx_regressors,:) = bootcoeff;
            boot_coeff{i} = boot_coeff{i}./S.factor_range(:);
            b(:,i) = b(:,i)./S.factor_range(:);
            cil(:,i) =  2*b(:,i) - quantile(boot_coeff{i}', 0.975)';
            ciu(:,i) =  2*b(:,i) - quantile(boot_coeff{i}', 0.025)';
        end
        S.T_mdl.b(ind_m,:) = median(b,2);
        S.T_mdl.cil(ind_m,:) = median(cil,2);
        S.T_mdl.ciu(ind_m,:) = median(ciu,2);
        boot.b = b;
        boot.cil = cil;
        boot.ciu = ciu;
        boot.coeff = boot_coeff;
        S.T_mdl.boot(ind_m) = boot;
    end
end
%% FIT_ED
% fit an exponential function using FMINCON
%
%=INPUT
%
%   Xtrain
%       The columns of the design matrix with values for the regressors,
%       for training.
%
%   Xttrain
%       The columns of the design that are the product of values of the
%       regressors multiplied to the numbers of days elapsed, for training.
%
%   ytrain
%       The observed unit count, for training.
%
%=OUTPUT
%   
%   b
%       coefficients
function b = fit_ed(X, Xt, t, y)
    opts=optimset('fmincon');
    opts.TolFun = 1e-10;
    opts.Display = 'off';
    % ones for the regressors associated with initial count (y0) and 0's for
    % the regressors associated with the decay rate (k)
    n_X = size(X,2);
    n_Xt = size(Xt,2);
    b0 = [ones(n_X,1); ones(n_Xt,1)]; 
    % A*b >= 0
    % The weighted sum of the decay terms has to be positive.
    if n_Xt > 1
        A_Xt = double(dec2bin(0:2^(n_Xt-1)-1)=='1');
        A_Xt = [ones(size(A_Xt,1),1), A_Xt];
    else
        A_Xt = 1;
    end
    ncol_A = size(A_Xt,1);
    A = [zeros(ncol_A,n_X), -A_Xt];
    b = fmincon(@(b) ls_exp(b, X, Xt, t, y), b0, A, zeros(ncol_A,1), ...
                            [], [], [], [], [], ...
                            opts);
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