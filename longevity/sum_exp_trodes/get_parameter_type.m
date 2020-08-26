%% GET_PARAMETER_TYPE
%   Extract the type of the model parameter
%=INPUT
%
%   parameter_name
%       A char vector that begins with either 'N1_' or 'k_'
%
%=OUTPUT
%   
%   parameter_type
%       A char vector
function parameter_type = get_parameter_type(parameter_name)
    if contains(parameter_name, 'N1f_')
        parameter_type = 'N1f';
    elseif contains(parameter_name, 'N1s_')
        parameter_type = 'N1s';
    elseif contains(parameter_name, 'N1_')
        parameter_type = 'N1';
    elseif contains(parameter_name, 'k_') && ~contains(parameter_name, 'N1_')
        parameter_type = 'k';
    elseif any(parameter_name==["kf", "ks"])
        parameter_type = 'k';
    else
        error('Cannot identify the type of parameter.')
    end
end