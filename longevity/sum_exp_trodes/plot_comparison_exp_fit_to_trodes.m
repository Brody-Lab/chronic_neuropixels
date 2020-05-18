% plot_comparison_exp_fit_to_trodes
from_scratch = false;
P = get_parameters;
add_folders_to_path
if from_scratch
    if ~exist('Cells', 'var')
        load(P.Cells_path);
    end
    S = compare_model_types_for_trodes(Cells, 'MCReps', 100)
    save([P.data_folder_path filesep 'comparison_exp_fit_to_trodes.met'], 'S')
else
    if ~exist('S', 'var')
        load([P.data_folder_path filesep 'comparison_exp_fit_to_trodes.met'])
    end
end

unique_fit_types = unique(S.fit_type);
xticklabel = cellfun(@(x) x(x~='_'), unique_fit_types, 'uni', 0);
n = numel(unique_fit_types);
figure
set(gca, P.axes_properties{:}, ...
         'XLim', [0.5, n+0.5], ...
         'XTick', 1:n, ...
         'XTickLabel', xticklabel)
for i = 1:n
    med = median(S.mse(S.fit_type==unique_fit_types(i)));
    ci=quantile(S.mse(S.fit_type==unique_fit_types(i)), [0.025, 0.975]);
    ci = 2*med-ci; % pivotal CI
    plot(i, med, 'ko', 'linewidth', 2)
    plot(i*[1,1], ci, 'k-', 'linewidth', 1)
end
ylabel('Mean square error')

% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'comparison_exp_fit_to_trodes'], P.figure_image_format{i})
end