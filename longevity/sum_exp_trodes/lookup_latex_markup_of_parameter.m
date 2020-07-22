function markup = lookup_latex_markup_of_parameter(parameter_name)
% LOOKUP_LATEX_MARKUP_OF_PARAMETER get the markup in Latex for axes labels
%
%=INPUT
%
%   parameter_name
%       A char vector
%
%=OUTPUT
%
%   markup
%       A char vector
    switch parameter_name
        case 'N1_const'
            markup = '$\beta_{0}^{N_{1}}$';
        case 'N1_AP_gt0'
            markup = '$\beta_{AP>0}^{N_{1}}$';
        case 'N1_DV_gtn2'
            markup = '$\beta_{DV>-2}^{N_{1}}$';
        case {'N1_AP', ...
              'N1_DV', ...
              'N1_ML', ...
              'N1_SO', ...
              'N1_SP', ...
              'N1_age', ...
              'N1_use'}
            markup = ['$\beta_{' get_regressor_name(parameter_name) '}^{N_{1}}$'];
        case 'k_AP_gt0'
            markup = '$\beta_{AP>0}^{k_{1}}$';
        case 'k_DV_gtn2'
            markup = '$\beta_{DV>-2}^{k_{1}}$';
        case {'k_AP', ...
              'k_DV', ...
              'k_ML', ...
              'k_SO', ...
              'k_SP', ...
              'k_age', ...
              'k_use'}
            markup = ['$\beta_{' get_regressor_name(parameter_name) '}^k$'];
        case 'alpha'
            markup = '$\alpha$';
        case 'k_slow'
            markup = '$\beta_{slow}$';
        case 'k_fast'
            markup = '$\beta_{fast}$';
        otherwise
            error('%s is an unknown parameter', char(parameter_name))
    end