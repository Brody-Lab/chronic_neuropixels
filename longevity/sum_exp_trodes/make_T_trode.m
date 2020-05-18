% MAKE_T_TRODE A table whose rows corresponds to electrodes from all
% recordings
%
%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILES
%
%=OPTIONAL INPUT
%
%   exclude_holderless
%       A logical scalar specifying whether exclude animals implanted
%       without a holder
%
%   x0
%       The first day after surgery to be examined.
%
%=OUTPUT
%
%   T_trode
%       A table whose rows corresponds to electrodes from all
%       recordings
%
%   factor_range
%       The range of the experimental factors
%
%   factor minima
%       The minima of the experimental factors
function [T_trode, exp_factors] = make_T_trode(Cells, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'exclude_holderless', true, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'mean_subtract_factors', false, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'normalize_factors', false, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'x0', P.x0, @(x) isscalar(x) && isnumeric(x))
addParameter(parseobj, 'unit_distance', P.unit_distance, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
T_trode = struct;
k = 0;
for i =1:numel(Cells)
    if Cells{i}.days_since_surgery < P_in.x0
        continue
    end
     if P_in.exclude_holderless && ...
       (Cells{i}.rat=="T170"||Cells{i}.rat=="T173")
        continue
    end
    k = k + 1;
    n_trode = sum(Cells{i}.electrodes.in_brain);
    T_trode.days_elapsed{k,1} = repmat(Cells{i}.days_since_surgery, n_trode, 1);
    T_trode.AP{k,1} = Cells{i}.electrodes.AP;
    T_trode.AP_gt0{k,1} = Cells{i}.electrodes.AP > 0;
    T_trode.DV{k,1} = Cells{i}.electrodes.DV;
    T_trode.DV_gtn2{k,1} = Cells{i}.electrodes.DV > -2;
    T_trode.ML{k,1} = Cells{i}.electrodes.ML;
    T_trode.SP{k,1} = ceil(Cells{i}.electrodes.index/2)*2/100; % um
    T_trode.SO{k,1} = ones(n_trode,1)*(Cells{i}.shank_plane=="coronal");
    for j = 1:n_trode
        is_nearby = abs(Cells{i}.electrode - Cells{i}.electrodes.index(j)) <=P_in.unit_distance;
        T_trode.unit{k,1}(j,1) = sum(is_nearby);
        T_trode.single_unit{k,1}(j,1) = sum(is_nearby & ...
                                            Cells{i}.ks_good(:));
    end
    T_trode.rat{k,1} = repmat(string(Cells{i}.rat), n_trode, 1);
    T_trode.bank{k,1} = Cells{i}.electrodes.bank;
    T_trode.Cells_index{k,1} = repmat(i, n_trode, 1);
    T_trode.identifier{k,1} = repmat(string(Cells{i}.identifier), n_trode,1);
end
T_trode = structfun(@(x) vertcat(x{:}), T_trode, 'uni', 0);
T_trode.days_since_init = T_trode.days_elapsed - P_in.x0;
T_trode = struct2table(T_trode);

exp_factors.name = P.sum_exp_trodes.exp_factors;
exp_factors.minima = min(T_trode{:, exp_factors.name});
exp_factors.range = max(T_trode{:, exp_factors.name}) - ...
                 min(T_trode{:, exp_factors.name});
if P_in.normalize_factors
    T_trode{:, exp_factors.name} = T_trode{:, exp_factors.name} - ...
                                      min(T_trode{:, exp_factors.name});
    T_trode{:, exp_factors.name} = T_trode{:, exp_factors.name} ./ ...
                                    max(T_trode{:, exp_factors.name});
end