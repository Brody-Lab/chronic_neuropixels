% PLOT_BEHAVIORAL_PERFORMANCE Compare the behavioral performance of animals
% when they are training and when they are tethered for recording
function [] = plot_behavioral_performance

plot_psychometrics('tethered', 1);
plot_trial_count_comparison;
plot_psychometric_comparison('sens')
plot_psychometric_comparison('abs_bias')
plot_psychometric_comparison('lapse')
plot_psychometric_comparison('prct_correct')

