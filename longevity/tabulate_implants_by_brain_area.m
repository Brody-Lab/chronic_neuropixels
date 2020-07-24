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
addParameter(parseobj, 'exclude_holderless', true, @(x) isscalar(x) && islogical(x))
addParameter(parseobj, 'group_by_brain_area', true, @(x) isscalar(x) && islogical(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
if ~isempty(P_in.group_by_brain_area)
    Cells = standardize_brain_area_names(Cells);
end
T_gain_noise = readtable(P.gain_noise_log_path);
T_implant = readtable(P.implant_log_path);
T_elec = struct;
for i = 1:numel(Cells)
   n_elec = numel(Cells{i}.electrodes.index);
   rat_name = unique(string(Cells{i}.rat));
   probe_serial = str2double(Cells{i}.probe_serial);
   T_elec.rat{i,1} = repmat(rat_name, n_elec, 1);
   T_elec.probe_serial{i,1} = repmat(probe_serial, n_elec, 1);
   T_elec.shank_plane{i,1} = repmat(string(Cells{i}.shank_plane), n_elec, 1);
   T_elec.sess_date{i,1} = repmat(Cells{i}.sess_date, n_elec, 1);
   if isfield(Cells{i}, 'penetration')
       T_elec.mm_inserted{i,1} = repmat(Cells{i}.penetration.depth_inserted, n_elec,1);
   else
       idx = T_implant.rat==rat_name & ...
             T_implant.neuropixels_sn == probe_serial;
       if sum(idx)~=1
           error('Cannot identify unique implant details')
       end
       T_elec.mm_inserted{i,1} = repmat(T_implant.depth_mm(idx), n_elec, 1);
   end
   T_elec.age{i,1} = repmat(get_age_at_implant(rat_name, probe_serial), n_elec,1);
   T_elec.AP{i,1} = Cells{i}.electrodes.AP;
   T_elec.DV{i,1} = Cells{i}.electrodes.DV;
   T_elec.ML{i,1} = Cells{i}.electrodes.ML;
   used_previously = probe_serial == T_gain_noise.probe_SN & ...
                     Cells{i}.sess_date > T_gain_noise.date_explanted; 
   T_elec.n_prev_use{i,1} = repmat(sum(used_previously), n_elec,1);
   if P_in.group_by_brain_area
       T_elec.brain_area{i,1} = Cells{i}.electrodes.brain_area;
   else
       T_elec.brain_area{i,1}=strings(n_elec,1);
   end
end
T_elec = structfun(@(x) vertcat(x{:}), T_elec, 'uni', 0);
T_elec = struct2table(T_elec);
T_elec.rat=categorical(T_elec.rat);

if P_in.exclude_holderless
   T_elec = T_elec(T_elec.rat ~= 'T170' & T_elec.rat ~= 'T173', :);
end

[G,ID] = findgroups(T_elec(:, {'rat', 'probe_serial', 'brain_area'}));
T=struct;
T.rat=categorical(ID.rat);
T.probe_serial = ID.probe_serial;
for i = 1:max(G)
    T.shank_plane(i,1) = unique(T_elec.shank_plane(G==i));
    T.mm_inserted(i,1) = unique(T_elec.mm_inserted(G==i));
    T.n_elec(i,1) = round(sum(G==i)/numel(unique(T_elec.sess_date(G==i))));
    T.age(i,1) = unique(T_elec.age(G==i));
    T.AP(i,1) = mean(T_elec.AP(G==i));
    T.DV(i,1) = mean(T_elec.DV(G==i));
    T.ML(i,1) = mean(T_elec.ML(G==i));
    T.n_prev_use(i,1) = unique(T_elec.n_prev_use(G==i));
    T.brain_area(i,1) = unique(T_elec.brain_area(G==i));
end
T=struct2table(T);
T.shank_plane = categorical(T.shank_plane);
T.brain_area = categorical(T.brain_area);