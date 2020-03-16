function Log = PB_import_implant_log
% PB_IMPORT_IMPLANT_LOG
kPaths = PB_get_constant('path');
Log = readtable(kPaths.implant_log);
%% Formating
for the_fields = Log.Properties.VariableNames; fie = the_fields{:};
    if iscell(Log.(fie)) && all(cellfun(@isstr, Log.(fie)))
        Log.(fie) = string(Log.(fie));
    end
end
Log.implant_date = datetime(Log.implant_date, 'Format', 'yyyy-MM-dd');
Log.date_explanted = datetime(Log.date_explanted, 'Format', 'yyyy-MM-dd');
Log.acute = logical(Log.acute);
%% Sort by date
[~,I] = sort(Log.implant_date, 'descend');
Log = Log(I,:);
%% Compute DV_mm
% DV0 is the DV of the insertion site
Log.DV_mm = Log.DV0_mm + -Log.depth_mm .* cosd(Log.coronal_deg);
%% Compute ML_mm
% ML0 is the ML of the insertion site
no_ML = isnan(Log.ML_mm);
Log.ML_mm(no_ML) = Log.ML0_mm(no_ML) + Log.ML0_mm(no_ML).*sind(Log.coronal_deg(no_ML));
%% nomenclature
Log.hemisphere(Log.hemisphere == string('left')) = string('L');
Log.hemisphere(Log.hemisphere == string('right')) = string('R');