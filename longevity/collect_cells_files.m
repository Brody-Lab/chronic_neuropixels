% script which loads the cells files for Thomas' and Adrian's neuropixels recordings,
% adds a few relevant fields, and concatenates the data into a big cell
% array called "Cells" 

% assumes each recording is all bank 0 or all bank 1. Things may break if
% that's not true.

bank_electrodes = {[1:191,193:384],[385:575,577:768],[769:863,865:960]}; % list the electrodes corresponding to each bank

%% import adrian's cells files and run any code which is specific to these files
recordings_log_csv_file = 'recordings_log.csv'; % will need to be redownloaded if google sheet changes!
recordings_table = read_recordings_log(recordings_log_csv_file); % convert to matlab table
use = find(~ismissing(recordings_table.cells_file_figure2) & recordings_table.usedForFigure2==1); % find recordings with a cells file path specified and load those files
for i=1:length(use)
    tmp = load(recordings_table.cells_file_figure2(use(i)));
    fields = fieldnames(tmp);
    for f=1:length(fields)
        Cells_AGB(i).(fields{f}) = tmp.(fields{f});
    end
    Cells_AGB(i) = import_penetration(Cells_AGB(i));
    Cells_AGB(i).days_since_surgery = days(datetime(Cells_AGB(i).sess_date) - datetime(Cells_AGB(i).penetration.date_implanted));  
    Cells_AGB(i).unique_bank = unique(Cells_AGB(i).bank);    
    Cells_AGB(i).bank_electrode_depths = (Cells_AGB(i).penetration.depth_inserted-0.2) - bank_electrodes{Cells_AGB(i).unique_bank+1}./100;    
    regions=arrayfun(@(x)x.name(1),Cells_AGB(i).penetration.regions);
    Cells_AGB(i).region_names=cell(size(Cells_AGB(i).regions));
    valid_region = Cells_AGB(i).regions>0;
    Cells_AGB(i).region_names(valid_region) = regions(Cells_AGB(i).regions(valid_region));    
    fprintf('.');    
end

%% import thomas' cells file and run any code which is specific to these files
P = get_parameters;
T = readtable(P.Thomass_recordings_path);
cells_paths=T.Cells_path;
for i=1:length(cells_paths)
    tmp = load(cells_paths{i});
    fields = fieldnames(tmp);
    for f=1:length(fields)
        Cells_TZL(i).(fields{f}) = tmp.(fields{f});
    end
    fprintf('.');
    Cells_TZL(i).DV = -Cells_TZL(i).anatom_loc_mm.DV; % match adrian's sign convention
    Cells_TZL(i).region_names = cellfun(@char,num2cell(Cells_TZL(i).cell_area),'UniformOutput',false);
    if length(Cells_TZL(i).region_names)~=length(Cells_TZL(i).cell_area)
        error('');
    end
    for k=1:length(Cells_TZL(i).DV)
        Cells_TZL(i).ks_good(k) = is_ks_good(Cells_TZL(i).raw_spike_time_s{k}); % for each unit, calculate whether it's a good SU based on Kilosort's metric for refractoriness
    end   
    Cells_TZL(i).probe_serial = num2str(Cells_TZL(i).meta.ap_meta.imDatPrb_sn )   ;
    Cells_TZL(i).rat=Cells_TZL(i).rat(1,:);        
    Cells_TZL(i).days_since_surgery = days(datetime(Cells_TZL(i).sess_date(1,:),'InputFormat','yyyy_MM_dd') - get_date_of_implant(Cells_TZL(i).rat,str2num(Cells_TZL(i).probe_serial)));     
    Cells_TZL(i).unique_bank = unique(Cells_TZL(i).bank);    
    implant_log = readtable(P.implant_log_path);
    row_idx = strcmp(implant_log.rat,Cells_TZL(i).rat) & implant_log.neuropixels_sn==str2num(Cells_TZL(i).probe_serial);        
    Cells_TZL(i).bank_electrode_depths = implant_log.depth_mm(row_idx) - bank_electrodes{Cells_TZL(i).unique_bank+1}./100;
    Cells_TZL(i).electrode = Cells_TZL(i).recording_site + 384*Cells_TZL(i).unique_bank;       
end

%% since Adrian's and Thomas' cells files have slightly different fields, concatenate into a cell array, not a struct array
Cells = [num2cell(Cells_AGB) num2cell(Cells_TZL)];


%% and then add any additional fields that don't require separate code for Thomas' and Adrian's sessions. 
% Thomas: This would be a good place to add new fields that seem useful for generating figures.
su_count=0;
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
    % Here I've decided to make a big struct array with one element per
    % cell, with fields specific to the cell and the parent session. But so
    % far I haven't found a use for such a struct in actually plotting
    % figures.
    Cells{i}.idx = i;    
    cell_fields = {'region_names','DV','electrode','bank','rat','days_since_surgery','probe_serial','idx'};
    good_cells = find(Cells{i}.ks_good(:) );
    for k=good_cells(:)'
        su_count = su_count+1;
        for c=1:length(cell_fields)
            if length(Cells{i}.(cell_fields{c}))>1 && ~ischar(Cells{i}.(cell_fields{c})) % element for each cell
                su_struct.(cell_fields{c})(su_count) = Cells{i}.(cell_fields{c})(k);
            else % element for each cells file
                if ischar(Cells{i}.(cell_fields{c}))
                    su_struct.(cell_fields{c}){su_count} = Cells{i}.(cell_fields{c});                
                else
                    su_struct.(cell_fields{c})(su_count) = Cells{i}.(cell_fields{c});                
                end
            end
        end
    end
end