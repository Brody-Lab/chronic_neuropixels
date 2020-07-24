function T = tabulate_implants_for_SU_curves(Cells)
% tabulate_implants_for_SU_curves
%
%   tabulate implants and sort the rows in a specificy way
P = get_parameters;
T=tabulate_implants_by_brain_area(Cells);

T(T.rat == 'T227' & T.probe_serial == 18194823381, :) = [];

T = T(~ismissing(T.brain_area),:);
T=T(ismember(string(T.brain_area), P.brain_areas), :);
T.area_idx = arrayfun(@(x) find(x==P.brain_areas), string(T.brain_area));
T.nAP = -T.AP;
T.nDV = -T.DV;
T=sortrows(T, {'area_idx', 'rat'});

% assign implant number as according to the row when an implant is first
% listed
[~,~,T.rat_idx] = unique(T.rat, 'stable');
[~,~,T.probe_serial_idx] = unique(T.probe_serial, 'stable');
T.implant = findgroups(T(:, {'rat_idx', 'probe_serial_idx'}));

T.AP = round(T.AP*10)/10;
T.DV = round(T.DV*10)/10;
T.ML = round(T.ML*10)/10;
T.new = T.n_prev_use > 0; 
T.n_prev_use = [];
T = T(:, {'brain_area', 'implant', 'rat', 'probe_serial', 'age', 'new', ...
          'mm_inserted', 'n_elec', 'shank_plane', 'AP', 'ML', 'DV'});

