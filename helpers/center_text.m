% CENTER_TEXT center align a text object horizontally or vertically in an
% axes
%
%=INPUT
%   text_obj
%       A text object
%
%   dim
%       A scalar specifying the axis of alignment.
%           1 - horizontal
%           2 - vertical
%           [1,2] - both

function [] = center_text(text_obj, dim)
assert(isa(text_obj, 'matlab.graphics.primitive.Text'));
ax = get(text_obj, 'Parent');
if any(dim==1)
    x_lim=get(ax, 'XLim');
    x0 = min(x_lim);
    x1= max(x_lim);
    if strcmpi(get(ax, 'XScale'), 'log')
        text_obj.Position(1) = (x1/x0/text_obj.Extent(3)*text_obj.Extent(1))^0.5;
    else
        text_obj.Position(1) = (x1-x0-text_obj.Extent(3)+text_obj.Extent(1))/0.5;
    end
elseif any(dim==2)
    y_lim=get(ax, 'YLim');
    y0 = min(y_lim);
    y1= max(y_lim);
    if strcmpi(get(ax, 'YScale'), 'log')
        text_obj.Position(2) = (y1/y0/text_obj.Extent(4)*text_obj.Extent(2))^0.5;
    else
        text_obj.Position(2) = (y1-y0-text_obj.Extent(4)+text_obj.Extent(2))/0.5;
    end
end