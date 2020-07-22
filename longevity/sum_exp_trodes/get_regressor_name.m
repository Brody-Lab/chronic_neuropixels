%% GET_REGRESSOR_NAME
%   Extract the name of the regressor from the name of the model parameter
%=INPUT
%
%   parameter_name
%       A char vector that begins with either 'N1_' or 'k_'
%
%=OUTPUT
%   
%   regressor_name
%       A char vector
function regressor_name = get_regressor_name(parameter_name)
    if get_parameter_type(parameter_name)=="N1"
        regressor_name = parameter_name(4:end);
    elseif get_parameter_type(parameter_name)=="k"
        regressor_name = parameter_name(3:end);
    else
        error('Cannot identify the type of parameter.')
    end
end