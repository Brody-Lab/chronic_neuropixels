function[]=figure5(varargin)
% FIGURE5 Make figure 6 from the manuscrip written by Luo, Bondy, et al.
%
% The figure shows that while tethered for Neuropixel recording without a
% cable commutator, rats perform a cognitively demanding task at a level
% similar to when they are untethered
%
%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from the
%       Brody Lab SQL repository.
parser = inputParser;
addParameter(parser, 'from_scratch', false, @(x) isscalar(x) && islogical(x))
parse(parser, varargin{:});
P_in = parser.Results;
    
add_folders_to_path
P = get_parameters;

if P_in.from_scratch
    sum_exp_trodes_run_mdl
else
    load(P.sum_exp_trodes.data_path);
end
sum_exp_trodes_plot_coef(S, 'mdl_inds', 1:5, 'rescale_parameters', false)