% GET_RANGE_OF_TRODE_FACTOR get the range of the values of experimental
% factors such as AP across electrodes
%
%=INPUT
%
%   T_trode
%       Table made using MAKE_T_TRODE
%
%=OUTPUT
%
%   factor_range
%       A row array whose length is NUMEL(P.ED_trode_regressors_all)
function factor_range = get_range_of_trode_factor(T_trode)
P=get_parameters;
idx_mm = ~contains(P.ED_trode_regressors_all, {'cat', 'SPA'}) & ...
          contains(P.ED_trode_regressors_all, {'_'});
f_range = range(T_trode{:,P.ED_trode_factors});

for i = 1:numel(P.ED_trode_regressors_all)
    if ~contains(P.ED_trode_regressors_all(i), {'_'}) || ...
       contains(P.ED_trode_regressors_all(i), {'cat', 'SPA'})
        factor_range(i) = 1;
    else
        name = P.ED_trode_regressors_all{i};
        name = strrep(name, 'y0_', '');
        name = strrep(name, 'k_', '');
        idx = strcmp(name,P.ED_trode_factors);
        if sum(idx) ~=1
            error('Cannot identify the correct factor name')
        end
        factor_range(i) = f_range(idx);
    end
end