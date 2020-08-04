function Cells = postprocess_Cells(Cells)
% POSTPROCESS_CELLS A script for calculating various metrics for Cells
% after collecting the Cells files
%% A230's AP coordnate got messed up
for i = 1:length(Cells)
    if Cells{i}.rat=="A230"
        Cells{i}.AP = 0.8 * ones(size(Cells{i}.AP));
    end
end
%% Standardize the datetime format
for i = 1:length(Cells)
    if isdatetime(Cells{i}.sess_date)
        continue
    else
        Cells{i}.sess_date = Cells{i}.sess_date(1,:);
        if contains(Cells{i}.sess_date,'-')
            Cells{i}.sess_date = datetime(Cells{i}.sess_date, 'input', 'yyyy-MM-dd');
        else
            Cells{i}.sess_date = datetime(Cells{i}.sess_date, 'input', 'yyyy_MM_dd');
        end
    end
end
%%
for i=1:length(Cells)
    if Cells{i}.days_since_surgery<0
        error('Session date is before the surgery date. Something is wrong.');
    end
    Cells{i}.ks_good = logical(Cells{i}.ks_good);
    Cells{i}.identifier = [Cells{i}.rat,num2str(Cells{i}.probe_serial),num2str(Cells{i}.unique_bank)]; % a unique configuration to track over time is defined by having the same rat, probe and bank
    Cells{i}.n_units = length(Cells{i}.ks_good);
    Cells{i}.n_good_units = sum(Cells{i}.ks_good);
    Cells{i}.mean_DV = mean(Cells{i}.DV(Cells{i}.ks_good));   
    Cells{i}.median_electrode_depth = median(Cells{i}.bank_electrode_depths(Cells{i}.bank_electrode_depths>0));
    Cells{i}.n_electrodes_in_brain = sum(Cells{i}.bank_electrode_depths>0);
    Cells{i}.n_units_superficial = length(Cells{i}.ks_good(Cells{i}.DV<Cells{i}.median_electrode_depth));
    Cells{i}.n_good_units_superficial = sum(Cells{i}.ks_good(Cells{i}.DV<Cells{i}.median_electrode_depth));    
    Cells{i}.n_units_deep = length(Cells{i}.ks_good(Cells{i}.DV>Cells{i}.median_electrode_depth));
    Cells{i}.n_good_units_deep = sum(Cells{i}.ks_good(Cells{i}.DV>Cells{i}.median_electrode_depth));
    %n_spikes = sum(cellfun(@numel,Cells{i}.raw_spike_time_s));
end
%% Fix Thomas's calculation of cell position
implant_log = readtable(P.implant_log_path);
implant_log.rat = string(implant_log.rat);
for i = 1:numel(Cells)
    if Cells{i}.rat(1) ~= 'A'
        idx = strcmp(implant_log.rat, Cells{i}.rat) & ...
                     implant_log.neuropixels_sn==str2num(Cells{i}.probe_serial);
        anatom_loc_mm=NP_get_cell_anatom_loc(implant_log(idx,:), Cells{i}.dist_from_tip_um);
        Cells{i}.AP = anatom_loc_mm.AP;
        Cells{i}.ML = anatom_loc_mm.ML;
        Cells{i}.DV = anatom_loc_mm.DV;
    end
end
%% Get electrodes positions (only the ones in the brain)
bank_electrodes = {[1:191,193:384],[385:575,577:768],[769:863,865:960]}; % list the electrodes corresponding to each bank
for i = 1:numel(Cells)
    if Cells{i}.rat(1) == 'A'
        electrodes_depth_mm = Cells{i}.bank_electrode_depths(Cells{i}.bank_electrode_depths>0);        
        [DV,ML,AP] = dual_angle_geometries(Cells{i}.penetration.angle.ML,...
                                           Cells{i}.penetration.angle.AP,...
                                           electrodes_depth_mm);
        Cells{i}.electrodes.DV = DV(:);
        Cells{i}.electrodes.ML = ML(:)+Cells{i}.penetration.craniotomy_ML;
        Cells{i}.electrodes.AP = AP(:)+Cells{i}.penetration.craniotomy_AP;
        Cells{i}.electrodes.index = [1:191, 193:384]' + 384*Cells{i}.unique_bank;
        trode_area = strings(960,1);
        for j = 1:numel(Cells{i}.penetration.regions)
            trode_area(Cells{i}.penetration.regions(j).electrodes) = Cells{i}.penetration.regions(j).name{1};
        end
        Cells{i}.electrodes.brain_area = trode_area(Cells{i}.electrodes.index);
        Cells{i}.electrodes.in_brain = Cells{i}.bank_electrode_depths(:)>0;
    else
        if any(strcmp(Cells{i}.rat, {'T170', 'T173', 'T176'}))
            tip_um = 137;
        else
            tip_um = 195;
        end
        elec_dist_from_tip_um = (bank_electrodes{Cells{i}.unique_bank+1}-1)*10 + tip_um;
        if min(elec_dist_from_tip_um) > min(Cells{i}.dist_from_tip_um) || ...
           max(elec_dist_from_tip_um) < max(Cells{i}.dist_from_tip_um)
            error('mismatch')
        end
        in_brain = Cells{i}.bank_electrode_depths>0;
        idx = strcmp(implant_log.rat, Cells{i}.rat) & ...
                     implant_log.neuropixels_sn==str2num(Cells{i}.probe_serial);
        Cells{i}.electrodes=NP_get_cell_anatom_loc(implant_log(idx,:), elec_dist_from_tip_um(in_brain));
        Cells{i}.electrodes.in_brain = in_brain(:);
        Cells{i}.electrodes.index = [1:191, 193:384]' + 384*Cells{i}.unique_bank;
        Cells{i}.electrodes.brain_area = NP_get_region_of_electrode(implant_log(idx,:), ...
                                                                   'sites', [1:191, 193:384]', ...
                                                                   'bank', Cells{i}.unique_bank);
    end
    Cells{i}.electrodes.bank = ones(Cells{i}.n_electrodes_in_brain,1)*Cells{i}.unique_bank;
    Cells{i}.electrodes.index=Cells{i}.electrodes.index(Cells{i}.electrodes.in_brain);
    Cells{i}.electrodes.brain_area=Cells{i}.electrodes.brain_area(Cells{i}.electrodes.in_brain);
end
%% Standardize anatomical coordinates
for i = 1:numel(Cells)
    Cells{i}.DV = -abs(Cells{i}.DV);
    Cells{i}.ML =  abs(Cells{i}.ML);
    Cells{i}.electrodes.DV = -abs(Cells{i}.electrodes.DV);
    Cells{i}.electrodes.ML =  abs(Cells{i}.electrodes.ML);
end
%% Check
for i = 1:numel(Cells)
for c = {'AP', 'ML', 'DV'}
    coord=c{:};
    if min(Cells{i}.electrodes.(coord)) - min(Cells{i}.(coord))>eps || ...
       max(Cells{i}.(coord)) - max(Cells{i}.electrodes.(coord)) >0.001
            error('mismatch')
    end
end
end
%% firing rate
for i = 1:numel(Cells)
    if isfield(Cells{i}, 'rec')
        file_time_s = Cells{i}.rec.ap_meta.fileTimeSecs;
    else
        file_time_s = Cells{i}.meta.ap_meta.fileTimeSecs;
    end
    Cells{i}.fr = cellfun(@(x) numel(x)/file_time_s, Cells{i}.raw_spike_time_s);
end
%% Standardize region names format
for i = 1:numel(Cells)
    region_names = Cells{i}.region_names(:);
    region_names(cellfun(@isempty, region_names)) = {''};
    Cells{i}.region_names = string(region_names);
end
%% probe shank plane orientation
for i = 1:numel(Cells)
    switch Cells{i}.rat
        case {'A230', 'A241', 'A242', 'A243', 'T181', 'T182'}
            Cells{i}.shank_plane = "sagittal";
        otherwise
            Cells{i}.shank_plane = "coronal";
    end
end
%% Standardize names of brain areas
