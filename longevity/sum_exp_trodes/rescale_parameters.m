function [b, cil, ciu] = rescale_parameters(T_res, T_mdl, T_regressor)
% RESCALE_PARAMETERS Rescale the parameter values to reflect the original
% range of the regressors
%
%=INPUT
%
%   T_res
%       A table containing the results of the model fitting
%
%   T_mdl
%       A table specifying whether parameters were allowed to vary and
%       which were fixed to zero
%
%   T_regressor
%       A table containing the range of each regressor
%
%=OUTPUT
%
%   b
%       A matrix of model parameter estimates. Its size is n_model_variant
%       by n_parameter_value
%
%   cil
%       A matrix of the lower confidence bounds of the parameter estimates,
%       computed by bootstrapping
%
%   ciu
%       A matrix of the upper confidence bounds
%
    T_range = unstack(T_regressor(:,{'name', 'range'}), 'range', 'name');
    T_range.const = 1;
    regressor_range = [1,1, 1];  % the first three parameters are alpha, k_fast, k_slow
    for i = 1:numel(T_mdl.Properties.VariableNames)
        regressor_name = get_regressor_name( T_mdl.Properties.VariableNames{i});
        regressor_range = [regressor_range, T_range.(regressor_name)];
    end
    b = T_res.b_med ./ regressor_range;
    cil = T_res.cil ./ regressor_range;
    ciu = T_res.ciu ./ regressor_range;
end