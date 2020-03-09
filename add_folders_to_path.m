% ADD_FOLDERS_TO_PATH add the directories for all analyses required for the
% manuscript's analysis to MATLAB's path
%
% No input or output.

function[]=add_folders_to_path()
addpath(genpath(fileparts(mfilename('fullpath'))));