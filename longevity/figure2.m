% This is a wrapper for doing the postprocessing on the Cells files
from_scratch = false;
%% Get the data
if from_scratch
    collect_cell_files
    postprocess_Cells
end
%% Plot average across all conditions
close all
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    plot_average_longevity(Cells, 'metric', metric{:});
end
%% Plot average, conditioned on DV
close all
bin_edges = [-10, -2,  -1,  0];
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'condition_on', 'DV', ...
                                  'anatom_bin_edges', bin_edges, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3);
end
%% Plot average, conditioned on AP
close all
bin_edges = [-8, 0, 4];
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'condition_on', 'AP', ...
                                  'anatom_bin_edges', bin_edges, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3);
end
%% Plot average, conditioned on bank
close all
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'condition_on', 'bank', ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3);
end