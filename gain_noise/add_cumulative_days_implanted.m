% ADD_CUMULATIVE_DAYS_IMPLANTED add a column to the GAIN_NOISE_LOG.CSV the
% cumulative number of days implanted
function [] = add_cumulative_days_implanted
P=get_parameters;
T= readtable(P.gain_noise_log_path);
T.cumul_days_implanted = days(T.date_explanted -T.date_implanted);
T.probe_sn = cellfun(@(x) str2double(x(12:22)), T.recording_id);
[~,idx] = sort(T.date_explanted);
T = T(idx,:);
unique_probes = unique(T.probe_sn);
for i = 1:numel(unique_probes)
    idx = T.probe_sn == unique_probes(i);
    T.days_implanted(idx) = cumsum(T.cumul_days_implanted(idx));
end
writetable(T, P.gain_noise_log_path);