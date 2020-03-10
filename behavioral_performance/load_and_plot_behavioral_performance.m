% LOAD_AND_PLOT_BEHAVIORAL_PERFORMANCE fetch the data from the database,
% process them, and then plot the results. Use PLOT_BEHAVIORAL_PERFORMANCE
% to skip the data fetching and processing steps
function[]=load_and_plot_behavioral_performance()
add_folders_to_path
make_recording_sessions_table
make_behavior_table
analyze_performance_by_rat
plot_behavioral_performance