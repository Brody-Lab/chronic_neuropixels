% FIGURE3 make all the plots of the figure from
% the pre-made CSV tables or from data being fetched online
%
%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from the
%       Brody Lab SQL repository.
function[]=figure3(varargin)
P_input = inputParser;
addParameter(P_input, 'from_scratch', false, @(x) isscalar(x) && islogical(x))
parse(P_input, varargin{:});
P_input = P_input.Results;
add_folders_to_path
if P_input.from_scratch
    make_recording_sessions_table % a table indicating the dates of recording sessions
    make_behavior_table % a table indicating dates of recording and also training sessions
    analyze_performance_by_rat
end
P = get_parameters;
figure('Pos', [100, 50, 900, 1200])
k = 0;
n_col = 3;
n_row = 2;
ax_hdl = [];

k=k+1;
ax_hdl(k) = subplot(n_col, n_row,k);
plot_behavioral_comparison('trials_done', 'axes', gca, 'samp_size_pos', [0.1, 0.9]);
label_hdl(k)=label_panel(gca, P.panel_labels(k+1), 'FontSize', P.panel_label_font_size);

k=k+1;
ax_hdl(k) = subplot(n_col, n_row,k);
plot_behavioral_comparison('prct_correct', 'axes', gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k+1), 'FontSize', P.panel_label_font_size);

k=k+1;
ax_hdl(k) = subplot(n_col, n_row,k);
plot_average_psychometrics('axes', gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k+1), 'FontSize', P.panel_label_font_size);
label_hdl(k).Position(1)=label_hdl(k-2).Position(1);

k=k+1;
ax_hdl(k) = subplot(n_col, n_row,k);
plot_behavioral_comparison('sens', 'axes', gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k+1), 'FontSize', P.panel_label_font_size);
label_hdl(k-1).Position(2)=label_hdl(k).Position(2);

k=k+1;
ax_hdl(k) = subplot(n_col, n_row,k);
plot_behavioral_comparison('abs_bias', 'axes', gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k+1), 'FontSize', P.panel_label_font_size);

k=k+1;
ax_hdl(k) = subplot(n_col, n_row,k);
plot_behavioral_comparison('lapse', 'axes', gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k+1), 'FontSize', P.panel_label_font_size);

for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure3'], P.figure_image_format{i})
end