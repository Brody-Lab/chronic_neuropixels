function SX = make_separate_design_matrices(T_trode, model_parameters)
% MAKE_SEPARATE_DESIGN_MATRICES make a separate design matrix for each term
% of the sum-of-exponentials regression model
%
%=INPUT
%
%   T_trode
%       A table of regressors for each electrode, made using MAKE_T_TRODE
%
%   model_parameters
%       A list of model parameters
%
%=OUTPUT
%
%   SX
%       A structure whose fields are the design matrices

SX = struct;
for i = 1:numel(model_parameters)
    parameter_type = get_parameter_type(model_parameters{i});
    regressor_name = get_regressor_name(model_parameters{i});
    SX.(parameter_type).(model_parameters{i}) = T_trode.(regressor_name);
end
SX = structfun(@struct2table, SX, 'uni', 0);
assert(numel(model_parameters) == sum(structfun(@(x) size(x,2), SX)))
