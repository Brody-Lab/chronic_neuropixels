function [] = get_data_for_comparing_yields_between_electrode_maps(varargin)
%% GET_DATA_FOR_COMPARING_YIELDS_BETWEEN_ELECTRODE_MAPS
%
%   How much is the yield higher if a sparse electrode map, skipping every
%   other electrode, were used, instead of an electrode map in which either
%   bank 0 or bank 1 were sampled densely?
%
%=OPTIONAL INPUT, POSITIONAL
%
%   1) reassemble
%       a scalar that is a logical, 0, or 1 indicating whether to remake
%       the dataset.  
P = get_parameters;
if nargin < 1
    reassemble = false;
else
    reassemble = varargin{1};
end
recording_ids = [   "T262_2020_09_03_bank0"
                    "T262_2020_09_03_bank1"
                    "T262_2020_09_03_banks_01_column"
                    "T262_2020_09_04_bank0"
                    "T262_2020_09_04_bank1"
                    "T262_2020_09_04_banks_01_column"
                    "T262_2020_09_08_bank1"
                    "T262_2020_09_08_banks01_column"
                    "T262_2020_09_09_bank0"
                    "T262_2020_09_09_bank1"
                    "T262_2020_09_09_banks01_column"];
kPath = PB_get_constant('path');
T = [];
for i = 1:numel(recording_ids)
    if reassemble
        NP_make_Cells(recording_ids{i}, 'align_to_behavior', false)
    end
    Cells_path = get_tzl_Cells_path(recording_ids{i});
    Cells = load(Cells_path);
    Tmp = struct;
    Tmp.sess_date = datetime(string(Cells.sess_date), 'format', 'yyyy_MM_dd');
    Tmp.recording_id = repmat(string(Cells.recording_id), numel(Tmp.sess_date), 1);
    Tmp.bank = Cells.bank;
    Tmp.electrode = Cells.electrode;
    Tmp.is_single = Cells.is_single;
    Tmp.brain_area = Cells.cell_area;
    Tmp.unitCount = Cells.unitCount;

    implant_date = get_date_of_implant(Cells.rat(1,:), Cells.meta.ap_meta.imDatPrb_sn);

    Tmp.days_elapsed = days(Tmp.sess_date - implant_date);

    if contains(Cells.meta.ap_meta.imRoFile, 'column')
        Tmp.electrode_map = repmat("column", numel(Tmp.sess_date), 1);
    elseif contains(Cells.meta.ap_meta.imRoFile, 'bank0')
        Tmp.electrode_map = repmat("bank0", numel(Tmp.sess_date), 1);
    elseif contains(Cells.meta.ap_meta.imRoFile, 'bank1')
        Tmp.electrode_map = repmat("bank1", numel(Tmp.sess_date), 1);
    else
        error('cannot determine electrode map')
    end
    T = [T; struct2table(Tmp)];
end
T.recording_id = categorical(string(T.recording_id));
T.brain_area = categorical(T.brain_area);
T.electrode_map = categorical(T.electrode_map);
save(P.electrode_map_comparison_data_path, 'T')
end