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
%   model_parameters
%       A char vector, string array, or cell array of char specifying the
%       model parameters. Each parameter must be a member of
%       P.possible_model_parameters
%
%   normalize_regressors
%       A scalar logical specifying whether to normalize the regressors
%
%   x0
%       The first day after surgery to be examined.
%
%   unit_distance
%       The distance, in electrode index, between an electrode and a unit's
%       central location for that unit to be counted for that electrode. A
%       UNIT_DISTANCE==0 indicates that each unit is counted once and for
%       the electrode where it had the largest spike waveform amplitude.
%
%=OUTPUT
%
%   T_trode
%       A table whose rows corresponds to electrodes from all
%       recordings
%
%   T_regressor
%       A table whose each row is a regressor, and the columns indicate the
%       pre-normalized minima and range
function [T_trode, T_regressor] = make_T_trode(Cells, varargin)
    parseobj = inputParser;
    P = get_parameters;
    addParameter(parseobj, 'exclude_holderless', true, @(x) isscalar(x) && islogical(x))
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.default_model_parameters, ...
        @(x) all(ismember(x, P.sum_exp_trodes.possible_model_parameters)))
    addParameter(parseobj, 'normalize_regressors', false, @(x) isscalar(x) && islogical(x))
    addParameter(parseobj, 'x0', P.x0, @(x) isscalar(x) && isnumeric(x))
    addParameter(parseobj, 'unit_distance', P.unit_distance, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    T_regressor.name = get_regressor_names(P_in.model_parameters);
    if ismember('age', T_regressor.name) % the animal's age at implantation
        T_delivery_date = readtable(P.delivery_date_path);
        T_delivery_date = unstack(T_delivery_date, 'delivery_date', 'rat_name');
    end
    if ismember('use', T_regressor.name) % the number of times the probe was previously used
        T_gain_noise = readtable(P.gain_noise_log_path);
    end
    Cells = standardize_brain_area_names(Cells);
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
        if ismember('age', T_regressor.name)
            age_days = days(Cells{i}.sess_date - T_delivery_date.(Cells{i}.rat)) - ...
                       Cells{i}.days_since_surgery;
            T_trode.age{k,1} = ones(n_trode,1)*age_days;
        end
        if ismember('AP', T_regressor.name)
            T_trode.AP{k,1} = Cells{i}.electrodes.AP;
        end
        if ismember('AP_gt0', T_regressor.name)
            T_trode.AP_gt0{k,1} = Cells{i}.electrodes.AP > 0;
        end
        if ismember('DV', T_regressor.name)
            T_trode.DV{k,1} = Cells{i}.electrodes.DV;
        end
        if ismember('DV_gtn2', T_regressor.name)
            T_trode.DV_gtn2{k,1} = Cells{i}.electrodes.DV > -2;
        end
        if ismember('ML', T_regressor.name)
            T_trode.ML{k,1} = Cells{i}.electrodes.ML;
        end
        if ismember('SO', T_regressor.name)
            T_trode.SO{k,1} = ones(n_trode,1)*(Cells{i}.shank_plane=="coronal");
        end
        if ismember('SP', T_regressor.name)
            T_trode.SP{k,1} = ceil(Cells{i}.electrodes.index/2)*2/100; % um
        end
        if ismember('use', T_regressor.name)
            used_previously = str2double(Cells{i}.probe_serial) == T_gain_noise.probe_SN & ...
                               Cells{i}.sess_date > T_gain_noise.date_explanted;          
            T_trode.use{k,1} = repmat(sum(used_previously),n_trode,1);
        end
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
        T_trode.brain_area{k,1} = Cells{i}.electrodes.brain_area;
    end
    T_trode = structfun(@(x) vertcat(x{:}), T_trode, 'uni', 0);
    T_trode.rat = categorical(T_trode.rat);
    T_trode.identifier = categorical(T_trode.identifier);
    T_trode.brain_area = categorical(T_trode.brain_area);
    T_trode.days_since_init = T_trode.days_elapsed - P_in.x0;
    T_trode = struct2table(T_trode);

    T_regressor.minimum = min(T_trode{:, T_regressor.name});
    T_regressor.range = max(T_trode{:, T_regressor.name}) - ...
                        min(T_trode{:, T_regressor.name});
    if P_in.normalize_regressors
        T_trode{:, T_regressor.name} = T_trode{:, T_regressor.name} - ...
                                       T_regressor.minimum;
        T_trode{:, T_regressor.name} = T_trode{:, T_regressor.name} ./ ...
                                        T_regressor.range;
    end
    T_regressor = struct2table(structfun(@transpose, T_regressor, 'uni',0));
end
%% GET_REGRESSOR_NAMES
%   extract the list of regressors from the list of model parameters
%
%=INPUT
%   
%   model_parameters
%       A char vector, string array, or cell array of char
function regressor_names = get_regressor_names(model_parameters)
    assert(ischar(model_parameters) || ...
           isstring(model_parameters) || ...
           iscellstr(model_parameters), ...
           '"model_parameters" must be a char vector, string array, or cell array of char')
    model_parameters = string(model_parameters);
    regressor_names = [];
    P = get_parameters;
    for i = 1:numel(P.sum_exp_trodes.possible_regressors)
        r = P.sum_exp_trodes.possible_regressors{i};
        if any(ismember({['N1_' r], ['k_' r]}, model_parameters))
            regressor_names = [regressor_names, string(r)];
        end
    end
end