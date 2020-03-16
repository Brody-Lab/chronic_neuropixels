% MAKE_FIGURE_BEHAVIORAL_PERFORMANCE make all the plots of the figure from
% the pre-made CSV tables or from data being fetched online
%
%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from the
%       Brody Lab SQL repository.
function[]=make_figure_behavioral_performance(varargin)
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
plot_psychometrics('tethered',1)
plot_behavioral_comparison('trials_done')
plot_behavioral_comparison('sens')
plot_behavioral_comparison('abs_bias')
plot_behavioral_comparison('lapse')
plot_behavioral_comparison('prct_correct')