function [] = tabulate_implants()
% TABULATE_IMPLANTS Create a table of all the implants and save it as an
% CSV file
    P = get_parameters;
    T = struct;
    k = 0;
    %% Adrian's implants
    penetrations = get_npx_penetrations();
    for i = 1:numel(penetrations)
        if ~ismember(penetrations(i).rat, P.rats)
            continue
        end
        k = k + 1;
        T.rat(k,1) = string(penetrations(i).rat);
        T.probe_serial(k,1) = str2double(penetrations(i).serial);
        T.AP(k,1) = penetrations(i).craniotomy_AP;
        T.ML(k,1) = penetrations(i).craniotomy_ML;
        if penetrations(i).hemisphere == "left"
            T.ML(k) = T.ML(k)*-1;
        end
        T.mm_inserted(k,1) = penetrations(i).depth_inserted;
        T.angle_coronal_deg(k,1) = penetrations(i).angle.ML;
        T.angle_sagittal_deg(k,1) = penetrations(i).angle.AP;
        T.shank_plane_angle(k,1) = abs(90-penetrations(i).probe_orientation);
        T.date_implanted(k,1) = datetime(penetrations(i).date_implanted, ...
                                        'InputFormat', 'yyyy-MM-dd');
        T.age(k,1) = get_age_at_time_of_implant(T.rat(k), T.date_implanted(k));
    end
    %% Thomas's implants
    T_implant = readtable(P.implant_log_path);
    idx = contains(T_implant.implant_type, 'Neuropixels') & ...
          ismember(T_implant.rat, P.rats);
    T_implant = T_implant(idx,:);    
    for i = 1:size(T_implant,1)
        k =k + 1;
        T.rat(k,1) = string(T_implant.rat(i));
        T.probe_serial(k,1) = T_implant.neuropixels_sn(i);
        T.AP(k,1) = T_implant.AP_mm(i);
        if T_implant.reference(i) == "IA0"
            T.AP(k,1) = T.AP(k,1) - 9;
        end        
        T.ML(k,1) = T_implant.ML_mm(i);
        T.mm_inserted(k,1) = T_implant.depth_mm(i);
        % negative means the tip is more to the
        % animal's left than the insertion site
        T.angle_coronal_deg(k,1) = T_implant.coronal_deg(i) * ...
                                   (1-2*(T_implant.hemisphere(i) == "left"));
        T.angle_sagittal_deg(k,1) = T_implant.sagittal_deg(i) * -1; 
        T.shank_plane_angle(k,1) = T_implant.neuropixels_shank_angle_to_coronal(i);
        T.date_implanted(k,1) = T_implant.implant_date(i);
        T.age(k,1) = get_age_at_time_of_implant(T.rat(k), T.date_implanted(k));
    end
    T = struct2table(T);
    [~,I] = unique(T(:,{'rat', 'probe_serial'}));
    T = T(I,:);
    T.n_prev_use = zeros(size(T,1),1);
    for i = 1:size(T,1)
        T.n_prev_use(i) = T.n_prev_use(i) + ...
                          sum(T.probe_serial(i) == T.probe_serial & ...
                              T.date_implanted(i) > T.date_implanted);
    end
    T = sortrows(T, {'AP', 'rat'}, 'descend');
    T.implant_number = (1:size(T,1))';
    T = T(:, [end, 1:end-1]);
    writetable(T, P.implants_path);
end
%% Get_age_at_time_of_implant
function age = get_age_at_time_of_implant(rat_name, date_implanted)
    delivery_date = bdata(['select deliverydate from ratinfo.rats ' ...
                           'where ratname regexp "{S}"'], char(rat_name));
    age = days(date_implanted- delivery_date);
end
