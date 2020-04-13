% FIT_EXPONENTIAL Fit an exponential or power law or a sum of either. The
% likelihood is returned for doing the likelihood ratio test. Assuming
% proportional noise.
%
%=INPUT
%
%   x
%       independent data, vector
%
%   y
%       dependent data, vector
%=OUTPUT
%
%   R
%       A structure
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   fit_type
%       A char array specifying the fit type
%
%   fit_type
%       A char array specifying the fit type
function R = fit_exponential(x,y, varargin)
parseobj = inputParser;
addParameter(parseobj, 'fit_type', 'exponential', @(x) ismember(x, {'exponential', 'power'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;