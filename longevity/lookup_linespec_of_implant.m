function [marker, colr] = lookup_linespec_of_implant(Cells, rat_name, probe_serial, T)
% LOOKUP_LINESPEC_OF_IMPLANT Each implant is assigned a combination of
% color and marker spec
%
%=INPUT
%
%   Cells
%       The cell array of structure made by "collect_cells_files" and "postprocess_cells"
%
%   rat_name
%       A char vector, string, or cell string that specifies that
%       name of a rat
%
%   probe_serial
%       The serial number of the probe
%
%=OPTIONAL INPUT, POSITIONAL
%   T
%       An implants table for specifying the color/marker order

validateattributes(Cells, {'cell'},{})
assert(ischar(rat_name) || isstring(rat_name) || iscellstr(rat_name), ...
    'RAT_NAME must be a char vector, string, or cell string.')
validateattributes(probe_serial, {'numeric'},{'scalar'})

if nargin < 4 && ~isempty(T)
    T = get_metrics_from_Cells(Cells);
    [~,ID]=findgroups(T(:,{'rat', 'probe_serial'}));
    ID = flip(ID);
    i = find(ID.rat==rat_name & ID.probe_serial == probe_serial);
else
    i = unique(T.implant(T.rat == rat_name & ...
                         T.probe_serial==probe_serial));
end
if numel(i)~=1
    error('Cannot find unique implant')
end

color_order = get(0, 'DefaultAxesColorOrder');
marker_order = {'o-', '*--', '^-.'};

ncolr = size(color_order,1);
nmark = numel(marker_order);

color_inds = repmat((1:ncolr), 1, nmark);
marker_inds = repmat((1:nmark), 1, ncolr);

% shuffle the combinations of color and marker
inds = rot90(reshape(1:((ncolr+1)*(nmark)), ncolr+1, nmark));
inds(inds>ncolr*nmark)=[];

s = rng;
rng('default')
ridx=randperm(ncolr*nmark, ncolr*nmark);
rng(s)

color_inds = color_inds(ridx);
marker_inds = marker_inds(ridx);

i = mod(i-1, ncolr*nmark)+1;
colr = color_order(color_inds(i),:);
marker = marker_order{marker_inds(i)};
