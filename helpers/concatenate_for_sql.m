% CONCATENATE_FOR_SQL format a DATETIME or NUMERIC array for querying using
% SQL commands
%
%=INPUT
%   data_input
%       A DATETIME or numeric array
%
%=OUTPUT
%   cat_car
%       A character array
%
%=EXAMPLE
%
% dt = datetime('today')-10:datetime('today');
% date_char=concatenate_for_sql(dt);
% sessid=bdata(['select sessid from sessions where ratname="T250" and sessiondate in (' date_char ...
%               ') order by sessiondate']);
function [cat_car]=concatenate_for_sql(data_input)
if isdatetime(data_input)
    data_input =string(datestr(data_input, 'yyyy-mm-dd'));
elseif isnumeric(data_input)
    data_input = num2cell(data_input);
    data_input = cellfun(@num2str, data_input, 'uni', 0);
else
    error('DATA_INPUT must be a DATETIME or NUMERIC array.')
end
cat_car='';
for i=1:numel(data_input)
    cat_car=[cat_car ',"'  data_input{i}  '"'];
end
cat_car = cat_car(2:end);