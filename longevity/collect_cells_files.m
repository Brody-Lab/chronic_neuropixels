function Cells = collect_cells_files(mode)

% function which loads the cells files for Thomas' and Adrian's neuropixels recordings,
% adds a few relevant fields, and concatenates the data into a big cell
% array called "Cells" 

% assumes each recording is all bank 0 or all bank 1. Things may break if
% that's not true.

if nargin
    mode=validatestring(mode,{'uncurated','curated','uncurated_matching_curated'},'collect_cells_files','mode',1);
else
    mode='uncurated';
end

switch mode
    case 'curated'
        cells_file_field='curated_cells_file';
    otherwise
        cells_file_field='cells_file';
end

bank_electrodes = {[1:191,193:384],[385:575,577:768],[769:863,865:960]}; % list the electrodes corresponding to each bank
P = get_parameters;
%% import adrian's cells files and run any code which is specific to these files
 % convert to matlab table, will need to be redownloaded if google sheet changes!
recordings_table = read_recordings_log(P.Adrians_recordings_path);
if mode=="uncurated_matching_curated"
    use = find(~ismissing(recordings_table.curated_cells_file) & recordings_table.used_in_chronic_npx_ms==1); % find recordings with a cells file path specified and load those files
else
    use = find(~ismissing(recordings_table.(cells_file_field)) & recordings_table.used_in_chronic_npx_ms==1); % find recordings with a cells file path specified and load those files    
end
for i=1:length(use)
    tmp = load(recordings_table.(cells_file_field)(use(i)));
    fields = fieldnames(tmp);
    for f=1:length(fields)
        Cells_AGB(i).(fields{f}) = tmp.(fields{f});
    end
    Cells_AGB(i) = import_penetration(Cells_AGB(i));
    Cells_AGB(i).days_since_surgery = days(datetime(Cells_AGB(i).sess_date) - datetime(strrep(Cells_AGB(i).penetration.date_implanted,' ','')));  
    Cells_AGB(i).unique_bank = unique(Cells_AGB(i).bank);    
    Cells_AGB(i).bank_electrode_depths = (Cells_AGB(i).penetration.depth_inserted-0.2) - bank_electrodes{Cells_AGB(i).unique_bank+1}./100;    
    regions=arrayfun(@(x)x.name(1),Cells_AGB(i).penetration.regions);
    Cells_AGB(i).region_names=cell(size(Cells_AGB(i).regions));
    valid_region = Cells_AGB(i).regions>0;
    Cells_AGB(i).region_names(valid_region) = regions(Cells_AGB(i).regions(valid_region));    
    fprintf('.');    
end
if mode=="uncurated"
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
else
    Cells=num2cell(Cells_AGB);
end
%% trim unnecessary fields
trim_fields = {'raw_spike_time_s','meanWfGlobalRaw','waveformSim','Trials','spike_time_s'};
for i=1:length(Cells)
    for f=1:length(trim_fields)
        if isfield(Cells{i},trim_fields{f})
            Cells{i} = rmfield(Cells{i},trim_fields{f});
        end
        if isfield(Cells{i},'waveform') && isfield(Cells{i}.waveform,trim_fields{f})
            Cells{i}.waveform = rmfield(Cells{i}.waveform,trim_fields{f});            
        end
    end
end
            