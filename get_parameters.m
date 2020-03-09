% GET_PARAMETERS Returns a structure with the parameters for analyses for
% the manuscript, such as file paths
%
%=OUTPUT
%
%   P
%       A structure with parameters for analyses
%
function P = get_parameters()
P.repository_path = fileparts(mfilename('fullpath'));
P.rat_info_path = [P.repository_path filesep 'behavioral_performance' filesep 'rat_info.csv'];
P.behavior_table_path = [P.repository_path filesep 'behavioral_performance' filesep 'behavior_table.csv'];
P.recording_sessions_path = [P.repository_path filesep 'behavioral_performance' filesep 'recording_sessions.csv'];
P.tzluo_path = [fileparts(P.repository_path) filesep 'tzluo'];