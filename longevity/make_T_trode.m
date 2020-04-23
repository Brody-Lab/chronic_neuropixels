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
function T_trode = make_T_trode(Cells, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'exclude_holderless', true, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x) && isnumeric(x))
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
    T_trode.DV{k,1} = Cells{i}.electrodes.DV;
    T_trode.ML{k,1} = Cells{i}.electrodes.ML;
    T_trode.SVP{k,1} = ceil(Cells{i}.electrodes.index/2)*2/100; % um
    T_trode.SPA{k,1} = ones(n_trode,1)*(Cells{i}.shank_plane=="coronal");
    for j = 1:n_trode
        is_nearby = abs(Cells{i}.electrode - Cells{i}.electrodes.index(j)) <=1; % two nearby sites
        T_trode.unit{k,1}(j,1) = sum(is_nearby);
        T_trode.single_unit{k,1}(j,1) = sum(is_nearby & ...
                                            Cells{i}.ks_good(:));
    end
    T_trode.rat{k,1} = repmat(string(Cells{i}.rat), n_trode, 1);
    T_trode.bank{k,1} = Cells{i}.electrodes.bank;
    T_trode.Cells_index{k,1} = repmat(i, n_trode, 1);
end
T_trode = structfun(@(x) vertcat(x{:}), T_trode, 'uni', 0);
T_trode = struct2table(T_trode);