% LABEL_PANEL Add a text annotation to the upper left hand corner of an
% axes
%
%=INPUT
%   ax
%       An AXES object
%
%   label
%       A char, cell, or string array of the annotation
%
%=OPTIONAL INPUT
%
%   FontSize
%       A numeric scalar specifying the font size of the label. The default
%       is the font size of the axes.
%
%   x
%       The left/right position of the annotation relative to the position
%       and width of the AXES. X=0 (default) places the annotation at the
%       left of the axes, and X=1 places it at the right of the axes.
%
%   y
%       The height of the annotation relative to the position and height of
%       the AXES. Y=1 (default) places the annotation at the top of the
%       axes, and Y=0 places it at the bottom of the axes.
%
%=OPTIONAL OUTPUT
%
%   1) The handle of the annotation object.
function varargout = label_panel(ax, label, varargin)
validateattributes(ax, {'matlab.graphics.axis.Axes'}, {'scalar'})
validateattributes(label, {'char', 'string', 'cell'}, {})
parseobj = inputParser;
addParameter(parseobj, 'FontSize', get(ax, 'FontSize'), ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}))
addParameter(parseobj, 'x', -0.1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}))
addParameter(parseobj, 'y', 1.1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}))
parse(parseobj, varargin{:});
P = parseobj.Results;
axes(ax);
outer_pos = get(ax, 'OuterPos');
pos(1) = outer_pos(1)+P.x*outer_pos(3); 
pos(2) = outer_pos(2)+P.y*outer_pos(4); 
pos(3:4) = [0,0]; % no line
pos = min(max(pos, 0),1);
hdl = annotation('textbox', pos, 'String', label,'FontSize', P.FontSize);
if nargout > 0
    varargout{1} = hdl;
end
