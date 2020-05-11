% GET_METRICS_FROM_CELLS get the neuronal stability metrics separated for
% each recording sessions and condition. 
%
%=INPUT
%
%   Cells
%       A cel array of structures made by COLLECT_CELL_FILES and
%       POSTPROCESS_CELLS
%
%=OUTPUT
%
%   T
%       A table with the metrics including the number of units, the number
%       of single units, event rate, and Vpp; the condition identifiers;
%       days elapsed; and identifying information for the recording.
%
%=OPTIONAL INPUT
%
%   AP_bin_edges
%       If CONDITION_ON includes 'AP', then this increasing vector
%       specifies the bin edges for binning AP positions
%
%   brain_area
%       If CONDITION_ON includes 'brain_area', then this char/cell/string
%       array specifies the brain areas that are conditions. For example,
%       to group the recordings from mPFC ('PrL', 'MO') and distinguish
%       that from other brain areas, then specify:
%           'brain_area', {{'PrL', 'MO'}, {'other'}}
%
%   condition_on
%       A char or cell array specifying the different conditions used for
%       averaging. The options are {'bank', 'ML', 'AP', 'DV',
%       'brain_area'}. Multiple conditions are allowed.
%
%   DV_bin_edges
%       If CONDITION_ON includes 'DV', then this increasing vector
%       specifies the bin edges for binning DV positions
%
%   exclude_3A
%       Exclude recordings from 3A probes.
%
%   exclude_holderless
%       A logical scalar specifying whether exclude animals implanted
%       without a holder
%
%   EI_bin_edges
%       If CONDITION_ON includes 'electode_index', then this increasing
%       vector specifies the bin edges for binning electrode indices. An
%       index of 1 indicates the electrode closest to the tip of the probe
%       shank, and an index of 960 indicates the electrode farthest from
%       the tip.
%
%   ML_bin_edges
%       If CONDITION_ON includes 'ML', then this increasing vector
%       specifies the bin edges for binning ML positions
%
%   x0
%       The first day to include.
function T = get_metrics_from_Cells(Cells, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'AP_bin_edges', P.AP_bin_edges, ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'brain_area', P.brain_area_groups, @(x) iscell(x)||ischar(x)||isstring(x))
addParameter(parseobj, 'condition_on', '', @(x) any(ismember(x, {'', ...
                                                             'electrode_index', ...
                                                             'EI', ...
                                                             'ML', ...
                                                             'AP', ...
                                                             'DV', ...
                                                             'shank_plane', ...
                                                             'brain_area'})))
addParameter(parseobj, 'DV_bin_edges', P.DV_bin_edges, ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'EI_bin_edges', P.EI_bin_edges, ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'ML_bin_edges', P.ML_bin_edges, ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'exclude_holderless', true, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'exclude_3A', false, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x) && isnumeric(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
%% Create the bin edges for conditions that aren't used
if ~contains(P_in.condition_on, 'AP')
    AP_bin_edges = [-inf, inf];
else
    AP_bin_edges = P_in.AP_bin_edges;
end
if ~contains(P_in.condition_on, 'DV')
    DV_bin_edges = [-inf, inf];
else
    DV_bin_edges = P_in.DV_bin_edges;
end
if ~contains(P_in.condition_on, 'ML')
    ML_bin_edges = [0, inf];
else
    ML_bin_edges = P_in.ML_bin_edges;
end
if ~any(contains(P_in.condition_on, 'electrode_index')) && ...
   ~any(contains(P_in.condition_on, 'EI'))
    EI_bin_edges = [0, 960];
else
    EI_bin_edges = P_in.EI_bin_edges;
end
if ~contains(P_in.condition_on, 'brain_area')
    brain_area = ""; % can't be a string because numel("")==1 while numel('')==0
else
    brain_area = cellfun(@string, P_in.brain_area, 'uni', 0);
    brain_area = cellfun(@(x) x(:)', P_in.brain_area, 'uni', 0);
    brain_area=brain_area(:)';
end
if ~any(contains(P_in.condition_on, 'shank_plane'))
    shank_plane = "";
else
    shank_plane = ["coronal", "sagittal"];
end
%% Create the results table
k = 0;
for i = 1:numel(Cells)
    if P_in.exclude_holderless && ...
       (Cells{i}.rat=="T170"||Cells{i}.rat=="T173")
        continue
    end
    if P_in.exclude_3A && ...
       Cells{i}.rat=="T176"
        continue
    end
    if Cells{i}.days_since_surgery < P_in.x0
        continue
    end
    % bin cells for the i-th recording
    AP_bin_cells = discretize(Cells{i}.AP, AP_bin_edges);
    DV_bin_cells = discretize(Cells{i}.DV, DV_bin_edges);
    ML_bin_cells = discretize(Cells{i}.ML, ML_bin_edges);
    EI_bin_cells = discretize(Cells{i}.electrode, EI_bin_edges);
    BA_bin_cells = ones(size(Cells{i}.region_names));
    for j = 1:numel(brain_area)
        if isempty(brain_area{j})
            continue
        elseif brain_area{j}=="other"
            idx = true(size(Cells{i}.region_names));
            for m = 1:numel(brain_area) 
                idx = idx & ~any(Cells{i}.region_names(:)==brain_area{m}, 2);
            end
            BA_bin_cells(idx)=j;
        else
            idx = any(Cells{i}.region_names(:)==brain_area{j}, 2);
            BA_bin_cells(idx)=j;
        end
    end
    SP_bin_cells = ones(size(Cells{i}.region_names));
    % bin electrodes for the i-th recording
    AP_bin_trode = discretize(Cells{i}.electrodes.AP, AP_bin_edges);
    DV_bin_trode = discretize(Cells{i}.electrodes.DV, DV_bin_edges);
    ML_bin_trode = discretize(Cells{i}.electrodes.ML, ML_bin_edges);
    EI_bin_trode = discretize(Cells{i}.electrodes.index, EI_bin_edges);
    BA_bin_trode = ones(numel(Cells{i}.electrodes.brain_area),1);
    for j = 1:numel(brain_area)
        if isempty(brain_area{j})
            continue
        elseif brain_area{j}=="other"
            idx = true(size(Cells{i}.electrodes.brain_area));
            for m = 1:numel(brain_area) 
                idx = idx & ~any(Cells{i}.electrodes.brain_area(:)==brain_area{m}, 2);
            end
            BA_bin_trode(idx)=j;
        else
            idx = any(Cells{i}.electrodes.brain_area(:)==brain_area{j}, 2);
            BA_bin_trode(idx)=j;
        end
    end
    SP_bin_trode = ones(size(Cells{i}.electrodes.brain_area));
    if numel(shank_plane) > 1
        SP_bin_cells = SP_bin_cells + (Cells{i}.shank_plane=="coronal");
        SP_bin_trode = SP_bin_trode + (Cells{i}.shank_plane=="coronal");
    end
    % keep only the cells that are associated with electrodes in the brain
    idx_in_brain = ismember(Cells{i}.electrode, Cells{i}.electrodes.index);
    % each row is a condition x recording
    c = 0; 
    for i_AP = 1:numel(AP_bin_edges)-1
    for i_DV = 1:numel(DV_bin_edges)-1
    for i_ML = 1:numel(ML_bin_edges)-1
    for i_EI = 1:numel(EI_bin_edges)-1
    for i_BA = 1:numel(brain_area)
    for i_SP = 1:numel(shank_plane)
        c=c+1;
        idx_cells = AP_bin_cells==i_AP & ...
                    DV_bin_cells==i_DV & ...
                    ML_bin_cells==i_ML & ...
                    EI_bin_cells==i_EI & ...
                    BA_bin_cells==i_BA & ...
                    SP_bin_cells==i_SP & ...
                    idx_in_brain;
        idx_trode = AP_bin_trode==i_AP & ...
                    DV_bin_trode==i_DV & ...
                    ML_bin_trode==i_ML & ...
                    EI_bin_trode==i_EI & ...
                    BA_bin_trode==i_BA & ...
                    SP_bin_trode==i_SP;
        if sum(idx_trode)<1
            continue
        end
        k=k+1;
        T.rat(k,1) = string(Cells{i}.rat);
        T.identifier(k,1) = string(Cells{i}.identifier);
        T.days_elapsed(k,1) = Cells{i}.days_since_surgery;
        T.condition(k,1) = c;
        T.i_AP(k,1) = i_AP;
        T.i_DV(k,1) = i_DV;
        T.i_ML(k,1) = i_ML;
        T.i_EI(k,1) = i_EI;
        T.i_BA(k,1) = i_BA;
        T.i_SP(k,1) = i_SP;
        T.AP_edges(k,1:2) = AP_bin_edges(i_AP:i_AP+1);
        T.DV_edges(k,1:2) = DV_bin_edges(i_DV:i_DV+1);
        T.ML_edges(k,1:2) = ML_bin_edges(i_ML:i_ML+1);
        T.EI_edges(k,1:2) = EI_bin_edges(i_EI:i_EI+1);
        
        T.trode_AP_avg(k,1) = mean(Cells{i}.electrodes.AP(idx_trode));
        T.trode_DV_avg(k,1) = mean(Cells{i}.electrodes.DV(idx_trode));
        T.trode_ML_avg(k,1) = mean(Cells{i}.electrodes.ML(idx_trode));
        T.trode_EI_avg(k,1) = mean(Cells{i}.electrodes.index(idx_trode));
        if ~isempty(brain_area{i_BA})
            T.brain_area(k,1) = string(strjoin(brain_area{i_BA}));
        else
            T.brain_area(k,1)="";
        end
        T.unit(k,1) = sum(idx_cells);
        T.single_unit(k,1) = sum(Cells{i}.ks_good(idx_cells));
        T.event_rate(k,1) = nansum(Cells{i}.fr(idx_cells));
        T.Vpp(k,1) = nanmean(Cells{i}.unitVppRaw(idx_cells));
        T.n_elec(k,1) = sum(idx_trode);
        T.shank_plane(k,1) = string(Cells{i}.shank_plane);
        T.probe_serial{k,1} = Cells{i}.probe_serial;
        T.sess_date(k,1) = Cells{i}.sess_date;
    end
    end
    end
    end
    end
    end
end
T.frac_single = T.single_unit./T.unit;
T.condition = findgroups(T.condition);
T.days_bin = discretize(T.days_elapsed, P.longevity_time_bin_edges);
T = struct2table(T);