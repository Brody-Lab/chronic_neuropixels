% A script to manually add a datafile that was excluded during initial
% processing
bank_electrodes = {[1:191,193:384],[385:575,577:768],[769:863,865:960]}; % list the electrodes corresponding to each bank
implant_log = readtable(P.implant_log_path);
if ~exist('Cells', 'var')
    load([P.data_folder_path filesep 'Cells.mat'])
end
i = 169;
file_path = ['X:\RATTER\PhysData\NP_sorted\Thomas\T212\T212_2019_09_03\spikesort_2020_03_15_15_00_46_ks2jrc' ...
             '\T212_2019_09_03_Cells.mat'];
tmp = load(file_path);
fields = fieldnames(tmp);
for f=1:length(fields)
    Cells{i}.(fields{f}) = tmp.(fields{f});
end
Cells{i}.DV = -Cells{i}.anatom_loc_mm.DV; % match adrian's sign convention
Cells{i}.region_names = cellfun(@char,num2cell(Cells{i}.cell_area),'UniformOutput',false);
if length(Cells{i}.region_names)~=length(Cells{i}.cell_area)
    error('');
end
for k=1:length(Cells{i}.DV)
    Cells{i}.ks_good(k) = is_ks_good(Cells{i}.raw_spike_time_s{k}); % for each unit, calculate whether it's a good SU based on Kilosort's metric for refractoriness
end   
Cells{i}.probe_serial = num2str(Cells{i}.meta.ap_meta.imDatPrb_sn )   ;
Cells{i}.rat=Cells{i}.rat(1,:);        
Cells{i}.days_since_surgery = days(datetime(Cells{i}.sess_date(1,:),'InputFormat','yyyy_MM_dd') - get_date_of_implant(Cells{i}.rat,str2num(Cells{i}.probe_serial)));     
Cells{i}.unique_bank = unique(Cells{i}.bank);    
row_idx = strcmp(implant_log.rat,Cells{i}.rat) & implant_log.neuropixels_sn==str2num(Cells{i}.probe_serial);        
Cells{i}.bank_electrode_depths = implant_log.depth_mm(row_idx) - bank_electrodes{Cells{i}.unique_bank+1}./100;
Cells{i}.electrode = Cells{i}.recording_site + 384*Cells{i}.unique_bank;       
%%
postprocess_Cells
%%
save([P.data_folder_path filesep 'Cells.mat'], 'Cells', '-v7.3')