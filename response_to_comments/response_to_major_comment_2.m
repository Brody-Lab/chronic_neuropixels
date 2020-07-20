% RESPONSE_TO_MAJOR_COMMENT_2 
%
% Plot the number of single units over time recorded from T181 separately
% for units recorded in (secondary) motor cortex and (anterior) dorsal
% striatum.
%
% "2. The results in Figure 2, especially the averages in Figure 2H,K, indicate
% severe losses in unit yield over time for probes implanted posterior of
% bregma and electrode sites in the dorsal 2mm of the rat brain. However,
% Figure 2B shows at least one animal (open circles) for which high neuron
% yield was obtained in motor cortex and dorsal striatum for at least 4
% months. First, is this from 1 or 2 probes? Whether it is from 1 or 2
% probes, the stability of the recording over time from day 1 is much
% greater than the other animals, and much better than what is expected
% from Figure 2H. Is the preservation of units over time for this animal
% due to the stability of units in dorsal striatum (presumably mostly >2mm
% below the surface) or also motor cortex?"
%% Load data
fprintf('\nLoading the variabe CELLS...')
if ~exist('Cells', 'var')
    load([P.data_folder_path filesep 'Cells.mat'])
end
fprintf('\nLoaded the variabe CELLS...')
%% Group data
P=get_parameters;
T = get_metrics_from_Cells(Cells, 'condition_on', {'EI', 'brain_area'}, ...
                                  'brain_area', {{'M2'}, {'Str'}}, ...
                                  'EI_bin_edges', 1*384 + [1, 384]);
T = T(T.rat=='T181',:);
T = T(T.n_elec > 1,:);
%% Plot data
figure('Position', [50,50,500,300]);
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
set(gca, 'xlim', [2^-0.5, 2^9.1], ...
         'xtick', 2.^(0:9));
for brain_area = ["M2", "Str"]
    idx = T.brain_area == brain_area;
    plot(T.days_elapsed(idx), T.single_unit(idx), 'o-', 'linewidth', 1);
end
legend({'Motor ctx', 'Dorsal striatum'}, 'Location', 'Best')
ylim([0,1].*ylim)
xlabel('Days since implant')
ylabel('Single units')
saveas(gcf, [P.plots_folder_path filesep mfilename], 'svg')