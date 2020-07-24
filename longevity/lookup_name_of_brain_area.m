function str = lookup_name_of_brain_area(brain_area)
% LOOKUP_NAME_OF_BRAIN_AREA provide the full description of a brain area
%
%=INPUT
%
%   BRAIN_AREA
%       A char vector of the abbreviation of a brain area
%=OUTPUT
%
%   str
%       A char vector of the full description of the brain area

assert(ischar(brain_area)||isstring(brain_area)||iscellstr(brain_area), ...
        '"brain_area" must be a char vector, string, or cell string')
switch brain_area
    case 'dmFC'
        str = 'Dorsomedial frontal cortex';
    case 'vmFC'
        str = 'Ventromedial frontal cortex';
    case 'MCtx'
        str = 'Motor cortex';    
    case 'dmStr'
        str = 'Dorsomedial striatum';        
    case 'vmStr'
        str = 'Ventromedial striatum'; 
    case 'piriform'
        str = 'Piriform areas';
    case 'S1'
        str = 'Primary somatosensory cortex';
    case 'lStr'
        str = 'Lateral striatum';
    case 'pallidum'
        str = 'Pallidum';
    case 'amygdala'
        str = 'Amygdaloid complex';
    case 'striatum tail'
        str = 'Striatum tail';
    case 'dpCtx'
        str = 'Dorsal-posterior cortex';
    case 'postsubiculum'
        str = 'Postsubiculum';
    case 'SC'
        str = 'Superior colliculus';
    case 'nIC'
        str = 'Posterior tectum';
    otherwise
        error('%s does not yet have a full description', char(brain_area))
end