% GET_CONDITION_NAMES get the names of conditions from the data table
%
%=INPUT
%
%   T
%       A table made using GET_METRICS_FROM_CELLS
%
%=OUTPUT
%
%   condition_names
%       A cell whose elements are char arrays
function condition_names = get_condition_names(T)
n_cond = numel(unique(T.condition));
if n_cond == 1
    condition_names = {'Overall_average'};
else
    condition_names = strings(1, n_cond);
    for i = 1:n_cond
        txt = '';
        
        unique_cond = unique(T.AP_edges, 'rows');
        if size(unique_cond,1)>1
            j = unique(T.i_AP(T.condition==i));
            txt = [txt ' AP [' num2str(unique_cond(j,1)) ', ' num2str(unique_cond(j,2)) '] mm'];
        end
        
        unique_cond = unique(T.DV_edges, 'rows');
        if size(unique_cond,1)>1
            j = unique(T.i_DV(T.condition==i));
            txt = [txt ' DV [' num2str(unique_cond(j,1)) ', ' num2str(unique_cond(j,2)) '] mm'];
        end
        
        unique_cond = unique(T.ML_edges, 'rows');
        if size(unique_cond,1)>1
            j = unique(T.i_ML(T.condition==i));
            txt = [txt ' ML [' num2str(unique_cond(j,1)) ', ' num2str(unique_cond(j,2)) '] mm'];
        end
        
        unique_cond = unique(T.EI_edges, 'rows');
        if size(unique_cond,1)>1
            j = unique(T.i_EI(T.condition==i));
            txt = [txt ' electrode [' num2str(unique_cond(j,1)) ', ' num2str(unique_cond(j,2)) ']'];
        end
        
        unique_cond = unique(T.brain_area);
        if numel(unique_cond)>1
            j = unique(T.i_BA(T.condition==i));
            txt = [txt ' ' unique_cond{j}];
        end
        
        unique_cond = unique(T.shank_plane);
        if numel(unique(T.i_SP))>1
            j = unique(T.i_SP(T.condition==i));
            txt = [txt ' ' unique_cond{j}];
        end
        
        condition_names{i} = txt(2:end);
    end
end