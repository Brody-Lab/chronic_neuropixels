function[] = plot_gain_noise_example()
P=get_parameters;
data_file_path = [P.gain_noise_fldr_path filesep P.gain_noise_example '.csv'];
D = readtable(data_file_path);

yticklabel = cellfun(@num2str, num2cell(0:5:P.noise_threshold_uV), 'uni', 0);
yticklabel{end} = ['\geq' num2str(P.noise_threshold_uV)];
figure('Position', P.figure_position_gn)
set(gca, P.axes_properties_behavior{:})
set(gca, 'Xlim', [0 961], ...
         'xtick', [], ...
         'YLim', [0, P.noise_threshold_uV], ...
         'Ytick', 0:5:P.noise_threshold_uV, ...
         'YtickLabel', yticklabel)
ylabel('RMS noise (uV)')
D.noise_uV(D.noise_uV>P.noise_threshold_uV)=P.noise_threshold_uV;
plot(D.noise_uV, 'ko', 'MarkerSize', 3);
T = readtable(P.gain_noise_log_path);
idx =strcmp(T.recording_id, P.gain_noise_example);
title([num2str(T.days_implanted(idx)) ' days implanted'], 'fontweight', 'normal')
% hdl_portion_imp = area([0 T.electrodes_implanted(idx)], max(ylim)*[1,1], ...
%                        'FaceAlpha', 0.1, ...
%                        'EdgeColor', 'none', ...
%                        'FaceColor', zeros(1,3));
xlabel('deep <- electrodes -> superficial')
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.gain_noise_fldr_path filesep 'RMS_noise_example'], P.figure_image_format{i})
end