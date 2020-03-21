%% and then add any additional fields that don't require separate code for Thomas' and Adrian's sessions. 
% Thomas: This would be a good place to add new fields that seem useful for generating figures.
% su_count=0;
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
    n_spikes = sum(cellfun(@numel,Cells{i}.raw_spike_time_s));
    switch Cells{i}.rat
        case {'T170', 'T173', 'T176'}
            Cells{i}.phase = '3A';
            Cells{i}.explantable = false;
        case {'T181', 'T182'}
            Cells{i}.phase = '3B';
            Cells{i}.explantable = false;
        otherwise
            Cells{i}.phase = '3B';
            Cells{i}.explantable = true;
    end
end
%% Fix Thomas's calculation of cell position
implant_log = readtable(P.implant_log_path);
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
    else
        if Cells{i}.phase == "3A"
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
    end
    Cells{i}.electrodes.bank = ones(Cells{i}.n_electrodes_in_brain,1)*Cells{i}.unique_bank;
end
%% Normalize anatomical coordinates
for i = 1:numel(Cells)
    Cells{i}.electrodes.DV = -abs(Cells{i}.electrodes.DV);
    Cells{i}.electrodes.ML =  abs(Cells{i}.electrodes.ML);
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