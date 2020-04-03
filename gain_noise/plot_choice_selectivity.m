% plot choice selectivity. In the first panel, it would show selectivity
% during the first penetration. In the second panel, it would show
% selectivity during the second animal. In the third panel, it would show
% selectivity during a recording from the third animal. I will show the
% heat maps all isolated units, which would provide information on the
% total number of isolated units recorded, which shouldn't be that
% different. I think I will take the first good session from each implant.
% The very first panel would actually be a diagram of the probe in anterior
% medial frontal cortex. Another possibility is to plot the average choice
% absolute choice selectivity for the three implants overlaid on each other
% in a summary plot to show that there is really not much difference
% between the signals. If that were true. If that were not true, then it
% wouldn't make sense to plot the averages. Instead, it might make more
% sense to plot just the heatmaps. This would not diminish my results on
% admFC when I finally publish them because there is no comparison here
% with other brain structures. Also, no modulation by stimulus intensity is
% shown.
%
% A) Schematic of the brain with the probe location
% B) heat maps showing the choice selectivity of isolated units from each
% recording
% C) summary plot averaging across the selectivity of isolated units

%% Fetch data from bucket
% 
% specify a path in the parameters script regarding where to get the files.
% Get both the Cells and the Trials file. Let's put them in the data
% folder. Ultimately, we are going to need to store the data on bucket.
% Until then, let's access them locally. We will likely create a bucket
% folder with the semi-processed data just for this manuscript. That would be fine with
% me. 
%% Sort trials, calculate the auROC

%% Separate function: plot the heatmaps

%% Separate function: plot the averages
% function [] = plot_choice_selectivity(varargin)
% parseobj = inputParser;
% addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
% parse(parseobj, varargin{:});

function [] = compute_choice_selectivity_examples
P = get_parameters;
if ~isfolder(P.choice_sel_data_path)
    error('Cannot find folder path')
end
for i = 1:numel(P.choice_sel_file_names)
    Cells_file_path = find_related_file(P.choice_sel_data_path, ...
                                        [P.choice_sel_file_names{i} ...
                                        '_Cells.mat'], ...
                                        'ignore_multiple', true);
    Trials_file_path = find_related_file(P.choice_sel_data_path, ...
                                        [P.choice_sel_file_names{i} ...
                                        '.Trials.mat'], ...
                                        'ignore_multiple', true);
    Cells_file_path=find_most_recent(Cells_file_path);
    Trials_file_path=find_most_recent(Trials_file_path);
    Cells{i} = load(Cells_file_path);
    S = load(Trials_file_path);
    if isfield(S, 'Trials')
        Trials{i} = S.Trials;
    else
        Trials{i} = S;
    end
end
%%
for i = 1:numel(P.choice_sel_file_names)
    vl_trials = ~Trials{i}.violated & ...
                  Trials{i}.responded & ...
                   Trials{i}.trial_type == 'a' & ...
                   Cells{i}.recorded{1};
    ChoiceMod{i} = PB_compute_choice_mod(Cells{i}, Trials{i}, ...
                                                   'compute_pval', false, ...
                                                   'ref', P.choice_sel_reference_event, ...
                                                   'steps_s', P.choice_sel_steps_s, ...
                                                   'bin_size_s',  P.choice_sel_bin_s, ...
                                                   'vl_trials', vl_trials);
end
%%
save(P.choice_sel_mat_path,'ChoiceMod')
%%
figure('Pos', [100,100,1000,350])
for i = 1:numel(P.choice_sel_file_names)
    subplot(1,numel(P.choice_sel_file_names), i)
    PB_plot_choice_mod(ChoiceMod{i}.time_s, ChoiceMod{i}.AUC, ...
                       'ax', gca, ...
                       'color_map_range', [0 1], ...
                       'color_bar', false)
    plot([0,0],ylim, 'k-', 'linewidth', 0.5)
    xlabel('Time from movement')
end