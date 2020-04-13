% GET_DETAILS_OF_RECORDINGS_INVENTORY get information on the sorted data
% files
%
%=INPUT
%   
%   T
%      A table listing the recordings with minimal information
%
%   sorted_data_path
%       path to the folder that contains all the spikesorted data
%
%=OUTPUT
%
%   T
%       A table listing the recordings with more information
function T = get_details_of_recordings_inventory(T, sorted_data_path)

n_sessions = size(T,1);
T.is_processed = false(n_sessions,1);
T.has_Cells = false(n_sessions,1);
T.has_Trials = false(n_sessions,1);
T.sorted_using = strings(n_sessions,1);
T.sorted_path = strings(n_sessions,1);
T.n_units = nan(n_sessions,1);
T.n_units_each_area = strings(n_sessions,1);
list_data_fldrs = dir(sorted_data_path);
data_fldr_names = string({list_data_fldrs.name})';
T.probe_sn = nan(n_sessions,1);
for i =1:n_sessions
    experimenter = get_experimenter_of_rat(T.rat_name{i});
    sess_fldr_path = [sorted_data_path    filesep ...
                      experimenter          filesep ...
                      T.rat_name{i}         filesep ...
                      T.recording_id{i}];
    if ~isfolder(sess_fldr_path)
        continue
    end
    date_str = datestr(T.session_date(i), 'yyyy_mm_dd');
    jrc4_file = find_related_file(sess_fldr_path, 'ap_res.mat', ...
                                    'ignore_missing', true, ...
                                    'ignore_multiple', true, ...
                                    'search_subfolders', true);  
    if isempty(jrc4_file)
        continue
    end
    jrc4_file = find_most_recent(jrc4_file);
    sorted_fdlr_path = fileparts(jrc4_file); %"spikesport_yyyy_mm_dd_hh_MM_SS_ks2jrc
    T.sorted_using{i,1} = 'Kilosort2';
    T.sorted_path{i,1} = sorted_fdlr_path;
    Cells_file = find_related_file(sorted_fdlr_path, 'Cells.mat', ...
                                    'ignore_missing', true, ...
                                    'ignore_multiple', false, ...
                                    'search_subfolders', true);  
    T.Cells_path{i,1} = Cells_file;
    Trials_file = find_related_file(sorted_fdlr_path, 'Trials.mat', ...
                                    'ignore_missing', true, ...
                                    'ignore_multiple', false, ...
                                    'search_subfolders', true);
    T.Trials_path{i,1} = Trials_file;
    T.is_processed(i) = ~isempty(T.Cells_path{i}) && ~isempty(T.Trials_path{i});
    if ~isempty(T.Cells_path{i})
        fprintf('\n %s  is processed',  T.recording_id(i))
        fprintf('\n   Loading %s', T.Cells_path{i})
        Cells = load(Cells_file);
        T.n_units(i) = numel(Cells.raw_spike_time_s);
        % list the number of each area, in the format of 'area1(n1)
        % area(n2)
        areas_with_units = unique(Cells.cell_area);
        area_str = '';
        for j = 1:numel(areas_with_units)
            n_unit_this_area = sum(Cells.cell_area ==  areas_with_units{j});
            if j > 1
                area_str = [area_str, ' '];
            end
            area_str = [area_str, areas_with_units{j} '(' ...
                        num2str(n_unit_this_area) ')'];
        end
        T.n_units_each_area{i} = area_str;
        T.probe_sn(i,1) = unique(Cells.meta.ap_meta.imDatPrb_sn);
    end
end