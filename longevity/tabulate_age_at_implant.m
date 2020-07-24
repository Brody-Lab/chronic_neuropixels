function []=tabulate_age_at_implant()
% TABULATE_AGE_AT_IMPLANT Make a CSV table of the age of the animal at each
% implant
P = get_parameters;
load(P.Cells_path)
T_delivery_date = readtable(P.delivery_date_path);
T_delivery_date = unstack(T_delivery_date, 'delivery_date', 'rat_name');
rat_names = cellfun(@(x) x.rat, Cells, 'uni', 0);
probe_serials = cellfun(@(x) str2double(x.probe_serial), Cells);
T=struct;
[G,T.rat, T.probe_serial] = findgroups(rat_names(:), probe_serials(:));
T.rat=categorical(T.rat);
for i = 1:max(G)    
    idx = G ==i;
    age = cellfun(@(x) days(x.sess_date - T_delivery_date.(char(T.rat(i))) - ...
                            x.days_since_surgery), Cells(idx));
    T.age_day(i,1)=age(1);
end
T = struct2table(T);
writetable(T, P.age_at_implant_path)