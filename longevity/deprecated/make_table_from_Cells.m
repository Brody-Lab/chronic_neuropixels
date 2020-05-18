% MAKE_TABLE_FROM_CELLS Assemble a table with a subset of the fields of
% CELLS
%
%=INPUT
%
%   Cells
%       A cell array made by COLLECT_CELL_FILES.m
%
%=OPTIONAL INPUT
%
%   fields
%       A char, string, or cell array specifying that fields of Cells that
%       would be incorporated into the table
%
%=OUTPUT
%
%   T
%       A table
function T = make_table_from_Cells(Cells, varargin)
parseobj = inputParser;
addParameter(parseobj, 'fields', {'days_since_surgery';
                                  'probe_serial';
                                  'n_units';
                                  'n_good_units';
                                  'rat';
                                  'sess_date';
                                  'unique_bank'}, ...
    @(x) validateattributes(x, {'string', 'cell', 'char'}, {}))
parse(parseobj, varargin{:});
P = parseobj.Results;
P.fields = string(P.fields);
Cells = Cells(:);
for i = 1:numel(P.fields)
    try 
        T.(P.fields{i}) = cellfun(@(x) x.(P.fields{i}), Cells);
    catch
        T.(P.fields{i}) = cellfun(@(x) x.(P.fields{i}), Cells, 'uni', 0);
    end
    switch P.fields{i}
        case 'probe_serial'
            T.(P.fields{i}) = cellfun(@str2double, T.(P.fields{i}));
    end
end
T = struct2table(T);
T.Vpp = cellfun(@(x) mean(x.unitVppRaw), Cells);
T.fr = cellfun(@(x) sum(x.fr), Cells);