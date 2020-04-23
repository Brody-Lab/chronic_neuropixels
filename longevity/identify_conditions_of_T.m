% IDENTIFY_CONDITIONS_OF_T list the conditions in the data table made by
% GET_METRICS_FROM_CELLS
%
%=INPUT
%   
%   T
%       A table made by GET_METRICS_FROM_CELLS
%
%=OUTPUT
%   
%   cond_name
%       A string column vector
function cond_name = identify_conditions_of_T(T)
cond_name = [];
if size(unique(T.AP_edges, 'row'),1) > 1
    cond_name = [cond_name; "AP"];
end
if size(unique(T.DV_edges, 'row'),1) > 1
    cond_name = [cond_name; "DV"];
end
if size(unique(T.ML_edges, 'row'),1) > 1
    cond_name = [cond_name; "ML"];
end
if size(unique(T.EI_edges, 'row'),1) > 1
    cond_name = [cond_name; "EI"];
end
if numel(unique(T.brain_area)) > 1
    cond_name = [cond_name; "brain_area"];
end