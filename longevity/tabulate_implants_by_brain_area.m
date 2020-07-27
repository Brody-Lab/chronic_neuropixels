function T = tabulate_implants_by_brain_area(Cells, varargin)
% TABULATE_IMPLANTSBY_BRAIN_AREA
%
% the animal, age at time of implant, probe newness, probe tip depth,
% estimated number of electrodes recorded from in that region, and shank
% orientation; the (AP,ML,DV) coordinates of the approximate midpoint of the
% probe's location within that brain region
%
%=INPUT
%
%   Cells
%       A cell of single-element structure made using "collect_cells_file"
%       and "postprocess_cells"
%
%=OUTPUT
%
%   T
%       A table enumerating each implant and sorted by brain_area
%
%=OPTIONAL INPUT
%
%   exclude_holderless
%       Exclude probes without a holder
%
%   group_by_brain_area
%       Whether to group implants by brain area
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'group_by_brain_area', true, @(x) isscalar(x) && islogical(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
if ~isempty(P_in.group_by_brain_area)
    Cells = standardize_brain_area_names(Cells);
end
T_implants = readtable(P.implants_path);
T_elec = struct;
for i = 1:numel(Cells)
    rat_name = unique(string(Cells{i}.rat));
    probe_serial = str2double(Cells{i}.probe_serial);
    if rat_name == "T170" || rat_name =="T173" ...  no holder
       || (rat_name == "T227" && probe_serial==18194823381) % no recording
        continue
    end    
    idx = T_implants.rat== rat_name & ...
          T_implants.probe_serial == probe_serial;
    if sum(idx)~=1
        error('The implant in Cells{%i} is not documented in the implants table', i)
    end
    n_elec = numel(Cells{i}.electrodes.index);
    T_elec.rat{i,1} = repmat(rat_name, n_elec, 1);
    T_elec.probe_serial{i,1} = repmat(probe_serial, n_elec, 1);
    if T_implants.shank_plane_angle(idx) == 0
        assert(Cells{i}.shank_plane=="coronal")
    elseif T_implants.shank_plane_angle(idx) == 90
        assert(Cells{i}.shank_plane=="sagittal")
    end
    T_elec.shank_plane_angle{i,1} = repmat(T_implants.shank_plane_angle(idx), n_elec, 1);
    T_elec.sess_date{i,1} = repmat(Cells{i}.sess_date, n_elec, 1);
    T_elec.mm_inserted{i,1} = repmat(T_implants.mm_inserted(idx), n_elec,1);
    T_elec.age{i,1} = repmat(T_implants.age(idx), n_elec,1);
    T_elec.AP{i,1} = Cells{i}.electrodes.AP;
    T_elec.DV{i,1} = Cells{i}.electrodes.DV;
    T_elec.ML{i,1} = Cells{i}.electrodes.ML;
    T_elec.n_prev_use{i,1} = repmat(T_implants.n_prev_use(idx), n_elec,1);
    if P_in.group_by_brain_area
       T_elec.brain_area{i,1} = Cells{i}.electrodes.brain_area;
    else
       T_elec.brain_area{i,1}=strings(n_elec,1);
    end
    T_elec.implant_number{i,1}=repmat(T_implants.implant_number(idx), n_elec,1);
end
T_elec = structfun(@(x) vertcat(x{:}), T_elec, 'uni', 0);
T_elec = struct2table(T_elec);
T_elec.rat=categorical(T_elec.rat);
T_elec.brain_area=categorical(T_elec.brain_area);
T_elec = T_elec(ismember(T_elec.brain_area, categorical(P.brain_areas)), :);

[G,ID] = findgroups(T_elec(:, {'rat', 'probe_serial', 'brain_area'}));
T=struct;
T.rat=categorical(ID.rat);
T.probe_serial = ID.probe_serial;
for i = 1:max(G)
    T.implant_number(i,1) = unique(T_elec.implant_number(G==i));
    T.shank_plane_angle(i,1) = unique(T_elec.shank_plane_angle(G==i));
    T.mm_inserted(i,1) = unique(T_elec.mm_inserted(G==i));
    T.n_elec(i,1) = round(sum(G==i)/numel(unique(T_elec.sess_date(G==i))));
    T.age(i,1) = unique(T_elec.age(G==i));
    T.AP(i,1) = round(mean(T_elec.AP(G==i)),1);
    T.DV(i,1) = round(mean(T_elec.DV(G==i)),1);
    T.ML(i,1) = round(mean(T_elec.ML(G==i)),1);
    T.n_prev_use(i,1) = unique(T_elec.n_prev_use(G==i));
    T.brain_area(i,1) = unique(T_elec.brain_area(G==i));
end
T=struct2table(T);

T.area_idx = arrayfun(@(x) find(x==P.brain_areas), string(T.brain_area));
T = sortrows(T, {'area_idx', 'implant_number'});
T = T(:, {'brain_area', 'implant_number', 'rat', 'probe_serial', 'age', 'n_prev_use', ...
          'mm_inserted', 'n_elec', 'shank_plane_angle', 'AP', 'ML', 'DV'});
writetable(T, P.implants_by_area_path);