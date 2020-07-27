function [marker, colr] = lookup_linespec_of_implant(rat_name, probe_serial)
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

assert(ischar(rat_name) || isstring(rat_name) || iscellstr(rat_name), ...
    'RAT_NAME must be a char vector, string, or cell string.')
validateattributes(probe_serial, {'numeric'},{'scalar'})
P=get_parameters;
T = readtable(P.implants_path);
i = find(strcmp(T.rat,rat_name) & T.probe_serial == probe_serial);
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

color_inds = color_inds(inds);
marker_inds = marker_inds(inds);

i = mod(i-1, ncolr*nmark)+1;
colr = color_order(color_inds(i),:);
marker = marker_order{marker_inds(i)};
