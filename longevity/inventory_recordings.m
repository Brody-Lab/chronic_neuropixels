% INVENTORY_RECORDINGS Search through folders to look for Neuropixels or
% tetrode recordings.
%
%=INPUT
%
%   folder_paths
%       A char, cell, or string array
%
%=OUTPUT
%
%   T
%       A table
function T = inventory_recordings(folder_paths, varargin)
validateattributes(folder_paths, {'cell', 'string', 'char'},{})
input_parser = inputParser;
% addParameter(input_parser, 'output', '', @(x) all(ismember(x, {'char', 'datetime', 'string'})))
parse(input_parser, varargin{:});
P = input_parser.Results;
folder_paths = string(folder_paths);
T = struct;
rat_list = [];
for i = 1:numel(folder_paths)
   rat_list = [rat_list; dir(folder_paths{i})];
end
k = 0;
for i = 1:numel(rat_list)
    % is the folder named after a rat?
    rat_fldr_path = fullfile(rat_list(i).folder, rat_list(i).name);
    if isfolder(rat_fldr_path) && ...
       numel(rat_list(i).name) == 4 && ... 
       isstrprop(rat_list(i).name(1), 'alpha') && ...
       all(isstrprop(rat_list(i).name(2:4), 'digit'))
        sess_list = dir(rat_fldr_path);
        for j = 1:numel(sess_list)
            sess_fldr_path = fullfile(sess_list(j).folder, sess_list(j).name);
            if isfolder(sess_fldr_path) && ...
               ~strcmp(sess_list(j).name, '.') && ...
               ~strcmp(sess_list(j).name, '..')
                rec_file_list = dir([sess_fldr_path  filesep '**' filesep '*.*']); % recursive
                is_tetrodes = any(arrayfun(@(x) contains(x.name, '.ncs'), rec_file_list)) || ...
                               any(arrayfun(@(x) contains(x.name, '.ntt'), rec_file_list));
                is_neuropixels = any(arrayfun(@(x) contains(x.name, '.ap.bin'), rec_file_list));
                if ~is_tetrodes && ~is_neuropixels
                    warning('\nUnknown physiology files: %s', sess_fldr_path)
                end
                if is_tetrodes
                    fprintf('\n%s: tetrodes', sess_fldr_path);
                    sess_date = sess_list(j).name(1:10);
                    sess_date = datetime(sess_date, 'inputFormat', 'yyyy-MM-dd');                    
                    % get time duration
                    % get time duration
                    % 16 byte header, and each record consists of uint64,
                    % uint32, uint32, uint32, and 512*int16, which is
                    % equal to 8 + 4 + 4 + 4 + 2*512 = 1044
                    % see https://neuralynx.com/software/NeuralynxDataFileFormats.pdf
                    % https://neuralynx.com/news/techtips/neuralynx-file-formats
                    ncs_file_list =  dir([sess_fldr_path  filesep '**' filesep '*.ncs*']); % recursive3
                    if numel(ncs_file_list) > 0 % there might be just NTT files
                        n_records = (ncs_file_list(1).bytes - 16)/1044;
                        dur_min = n_records * 15728/1e6 / 60; 
                    else
                        dur_min = nan;
                    end
                else
                    fprintf('\n%s: neuropixels', sess_fldr_path);
                    if all(isstrprop(sess_list(j).name(6:9), 'digit')) % e.g. 2019
                        sess_date = sess_list(j).name(6:15);
                    else
                        sess_date = sess_list(j).name(12:21);
                    end
                    sess_date = datetime(sess_date, 'inputFormat', 'yyyy_MM_dd');
                    ap_file = dir([sess_fldr_path, filesep '**' filesep '*.ap.bin']); % recursive3  
                    % 384 * 1 byte x 2 + 1 byte digital word
                    dur_min = ap_file(1).bytes/770/3e4/60;
                end
                % log in the table
                k = k + 1;
                T(k,1).session_date = sess_date;
                T(k,1).rat_name = rat_list(i).name;
                T(k,1).recording_id = sess_list(j).name;
                T(k,1).recording_id = strrep(T(k,1).recording_id, '_g0', '');
                if is_tetrodes
                    T(k,1).sensor = 'tetrodes';
                else
                    T(k,1).sensor = 'Neuropixels';
                end
                T(k,1).duration_min = round(dur_min);
            end
        end
    end
end
T = struct2table(T);
T.rat_name = string(T.rat_name);
T.recording_id = string(T.recording_id);
T.sensor = string(T.sensor);