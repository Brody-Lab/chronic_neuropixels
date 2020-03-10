function[]=plot_trial_count_comparison(varargin)
P_input = inputParser;
addParameter(P_input, 'min_trials', 10, @(x) isscalar(x) && (x == 0||x==1))
parse(P_input, varargin{:});
P_input = P_input.Results;
P = get_parameters;
T = readtable(P.performance_by_session_path);
T = T(T.trials_done>P_input.min_trials,:); % when techs mis-click the run-session button
figure('Position', P.figure_position_behavioral_comparison)
set(gca, P.axes_properties_behavior{:}, ...
         'DataAspectRatio', [1,1,1], ...
         'XScale', 'log', ...
        'YScale', 'log', ...
        'XLim', [100, 2e3], ...
        'YLim', [100, 2e3])
unique_rats = unique(T.rat);
for i=1:numel(unique_rats)
    i_tethered = strcmp(T.rat, unique_rats{i}) & T.tethered;
    i_untether = strcmp(T.rat, unique_rats{i}) & ~T.tethered;
    data_tethered=T.trials_done(i_tethered);
    data_untether=T.trials_done(i_untether);
    med_tethered =nanmedian(data_tethered);
    med_untether = nanmedian(data_untether);
    plot(med_untether, med_tethered, 'o');
end
diag_coords = [min([xlim,ylim]), max([xlim,ylim])];
plot(diag_coords,diag_coords, 'k','linewidth', 0.5)
xlabel('Untethered')
ylabel('Tethered for recording')
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.behavior_plots_folder filesep 'trial count comparison'], P.figure_image_format{i})
end