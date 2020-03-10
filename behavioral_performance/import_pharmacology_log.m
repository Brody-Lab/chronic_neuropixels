function Log = PB_import_pharmacology_log()

warning('off', 'MATLAB:table:ModifiedVarnames')
P = get_parameters;
Log = readtable(P.pharmacology_log_path);
warning('on', 'MATLAB:table:ModifiedVarnames')

% extract dates of pharmacological manipulation
pharma_dates = datetime(Log.Date, 'Format', 'yyyyMMdd');
pharma_dates = pharma_dates(~isnat(pharma_dates)); % take only the real date-times
pharma_dates = str2double(string(pharma_dates));
pharma_dates = unique(pharma_dates);

% convert certain columns to string arrays to facilitate future manipulation
Log.Manipulation = string(Log.Manipulation);
Log.Hemisphere = string(Log.Hemisphere);
Log.Rat = string(Log.Rat);

Log.RightDoseMultipler(isnan(Log.RightDoseMultipler)) = 1;
Log.RightVolumeMultipler(isnan(Log.RightVolumeMultipler)) = 1;

% findthe correct brain_areas
Log.brain_area = repmat(string, size(Log.Rat,1),1);
rats_that_got_drugs = unique(Log.Rat);
Implant_log = PB_import_implant_log;
for i = 1:numel(rats_that_got_drugs)
    idx = Implant_log.rat == rats_that_got_drugs(i) & ...
          Implant_log.implant_type == string('cannula');
    if sum(idx)<1
        error('Cannot find implant information for %s', rats_that_got_drugs{i});
    end
    implanted_area = unique(Implant_log.brain_area(idx));
    if numel(implanted_area) > 1
        error('The code cannot handle cannulas in differet brain area')
    end
    idx_infused = Log.Rat == rats_that_got_drugs(i) & ...
                  (Log.Manipulation == 'muscimol' | ...
                  Log.Manipulation == 'saline'); % the other manipulation is 'isofluorane'
    Log.brain_area(idx_infused) = implanted_area;
end
Log.brain_area = string(Log.brain_area);