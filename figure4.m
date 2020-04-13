% FIGURE4 Make figure 4 from the manuscrip written by Luo, Bondy, et al.
%
%   The figure show that explanted probes and unimplanted probes have
%   similar input-referred noise and can acquire neural signals of similar
%   quality.

%% Load the cells data
from_scratch = false; % Do you want to reassemble the data files from scratch?
P=get_parameters;
if from_scratch
    collect_cell_files
    postprocess_Cells
else
    if ~exist('Cells', 'var')
        fprintf('Loading the variabe CELLS...')
        load([P.data_folder_path filesep 'Cells.mat'])
    end
end
%% Make the panels
figure4_ABCD
figure4_EFGH(Cells)
figure4_IJKL