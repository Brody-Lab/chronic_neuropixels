function [] = compare_yields_between_electrode_maps(varargin)
%% GET_DATA_FOR_COMPARING_YIELDS_BETWEEN_ELECTRODE_MAPS
%
%   How much is the yield higher if a sparse electrode map, skipping every
%   other electrode, were used, instead of an electrode map in which either
%   bank 0 or bank 1 were sampled densely?
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   reassemble
%       a scalar that is a logical, 0, or 1 indicating whether to remake
%       the dataset.  
%
%   SUs_only
%       Include only single units
parser = inputParser;
addParameter(parser, 'reassemble', false, @(x) isscalar(x) && (x==0 || x== 1))
addParameter(parser, 'SUs_only', false, @(x) isscalar(x) && (x==0 || x== 1))
parse(parser, varargin{:});
P_in = parser.Results;

P = get_parameters;
get_data_for_comparing_yields_between_electrode_maps(P_in.reassemble)
load(P.electrode_map_comparison_data_path, 'T');
if P_in.SUs_only
    T = T(logical(T.is_single),:);
end
C = groupcounts(T, {'days_elapsed', 'electrode_map'});
figure
set(gca, P.axes_properties{:}, ...
         'Xscale', 'log', ...
         'ActivePositionProperty', 'OuterPosition')
handle = [];
k = 0;
for emap = unique(C.electrode_map)'
    idx = C.electrode_map == emap;
    k = k +1;
    handle(k) = plot(C.days_elapsed(idx), C.GroupCount(idx), 'o-');
end
legend(handle, cellstr(unique(C.electrode_map)))
xlabel('Days since implant')
ylim(ylim.*[0,1])
if P_in.SUs_only
    ylabel('Single units')
else
    ylabel('Units')
end

for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'compare_yields_between_electrode_maps'], P.figure_image_format{i})
end