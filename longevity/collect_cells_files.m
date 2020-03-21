% script which loads the cells files for Thomas' and Adrian's neuropixels recordings,
% adds a few relevant fields, and concatenates the data into a big cell
% array called "Cells" 

% assumes each recording is all bank 0 or all bank 1. Things may break if
% that's not true.

bank_electrodes = {[1:191,193:384],[385:575,577:768],[769:863,865:960]}; % list the electrodes corresponding to each bank
P = get_parameters;
%% import adrian's cells files and run any code which is specific to these files
 % convert to matlab table, will need to be redownloaded if google sheet changes!
recordings_table = read_recordings_log(P.Adrians_recordings_path);
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
T = readtable(P.Thomass_recordings_path);
implant_log = readtable(P.implant_log_path);
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
    row_idx = strcmp(implant_log.rat,Cells_TZL(i).rat) & implant_log.neuropixels_sn==str2num(Cells_TZL(i).probe_serial);        
    Cells_TZL(i).bank_electrode_depths = implant_log.depth_mm(row_idx) - bank_electrodes{Cells_TZL(i).unique_bank+1}./100;
    Cells_TZL(i).electrode = Cells_TZL(i).recording_site + 384*Cells_TZL(i).unique_bank;       
end

%% since Adrian's and Thomas' cells files have slightly different fields, concatenate into a cell array, not a struct array
Cells = [num2cell(Cells_AGB) num2cell(Cells_TZL)];

%% 
