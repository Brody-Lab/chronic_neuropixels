function [] = list_Thomass_recordings

folder_paths = {'Y:\RATTER\PhysData\Raw\Thomas';
                'Z:\RATTER\PhysData\Raw\Thomas'};
sorted_data_path = 'X:\RATTER\PhysData\NP_sorted';
map_bucket_drive
map_archive_drive
map_archiveme_drive
T0 = inventory_recordings(folder_paths);
T1 = T0(T0.sensor == "Neuropixels" & ...
        T0.duration_min<30 & ...
        ~contains(T0.recording_id, {'anesthetized', 'explant'}), :);
T2 = get_details_of_recordings_inventory(T1, sorted_data_path);
T3 = T2(T2.sorted_path ~= "" & ...
        ~isempty(T2.Cells_path), :);
P = get_parameters;
if ~isfolder(P.longevity_folder_path)
        mkdir(P.longevity_folder_path)
end
writetable(T3, [P.longevity_folder_path filesep 'Thomass_recordings.csv'])