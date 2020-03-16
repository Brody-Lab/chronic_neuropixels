function[]=plot_behavioral_comparison(data_type)
P = get_parameters;
T = readtable(P.performance_by_rat_path);
figure('Position', P.figure_position_behavioral_comparison)
set(gca, P.axes_properties_behavior{:}, ...
         'DataAspectRatio', [1,1,1])
if isfield(P.custom_axes_properties, data_type)
    set(gca, P.custom_axes_properties.(data_type){:})
end
unique_rats = unique(T.rat);
for i=1:numel(unique_rats)
    i_tethered = strcmp(T.rat, unique_rats{i}) & T.tethered;
    i_untether = strcmp(T.rat, unique_rats{i}) & ~T.tethered;
    data_tethered(i)=T.(data_type)(i_tethered);
    data_untether(i)=T.(data_type)(i_untether);
    if strcmp(data_type, 'bias')
        data_tethered(i)=abs(data_tethered(i));
        data_untether(i)=abs(data_untether(i));
    end
    plot(data_untether(i), data_tethered(i), 'o', 'Color', P.color_order(i,:));
end
pval = num2str(signtest(data_tethered, data_untether));
diag_coords = [min([xlim,ylim]), max([xlim,ylim])];
plot(diag_coords,diag_coords, 'k','linewidth', 0.5)
title(['p = ' pval], 'FontSize', 10, 'FontWeight', 'normal')
xlabel('Untethered')
ylabel('Tethered for recording')
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.behavior_plots_folder filesep 'psychometric_' data_type], P.figure_image_format{i})
end