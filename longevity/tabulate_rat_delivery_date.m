function [] = tabulate_rat_delivery_date(Cells)
% TABULATE_RAT_DELIVERY_DATE
%   Create a CSV table of the deliver date of each rat
%
%=INPUT
%
%   Cells
%       A structure containing the physiology data created by
%       COLLECT_CELLS_FILES.m

    P=get_parameters;
    rat_name = unique(cellfun(@(x) x.rat, Cells, 'uni', 0));
    rat_name = cellfun(@(x) [x '|'], rat_name, 'uni', 0);
    rat_str = strcat(rat_name{:});
    rat_str = rat_str(1:end-1);
    T = struct;
    [T.rat_name, T.delivery_date] = ...
        bdata('select ratname, deliverydate from ratinfo.rats where ratname regexp "{S}"', rat_str);
    T=struct2table(T);
    writetable(T, P.delivery_date_path)
end