% a script for finding Thomas's data

sorted_data_path = {'X:\RATTER\PhysData\NP_sorted\Thomas';
                    'X:\RATTER\PhysData\NP_sorted\Chuck';};
T0 = inventory_recordings(sorted_data_path);
[~,I] = sort(T0.recording_id);
T0 = T0(I,:);   
T1 = T0(T0.sensor == "Neuropixels" & ...
        T0.duration_min<30 & ...
        ~contains(T0.recording_id, {'anesthetized', 'explant'}), :);
T2 = get_details_of_recordings_inventory(T1, 'X:\RATTER\PhysData\NP_sorted\');
T3 = T2(T2.sorted_path ~= "" & ...
        ~isempty(T2.Cells_path), :);
[~, I] = unique(T3.Cells_path);
T3 = T3(I,:);
[~,I] = sort(T3.recording_id);
T3 = T3(I,:);
P = get_parameters;
if ~isfolder(P.longevity_folder_path)
    mkdir(P.longevity_folder_path)
end
writetable(T3, P.Thomass_recordings_path)