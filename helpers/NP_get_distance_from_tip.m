function dist_um = NP_get_distance_from_tip(implant, varargin)

input_parser = inputParser;
addParameter(input_parser, 'bank', 0, @(x) (ischar(x) || isnumeric(x)))
addParameter(input_parser, 'sites', 1:384, @(x) isnumeric(x))    % electrode index
parse(input_parser, varargin{:});
bank = input_parser.Results.bank;
sites = input_parser.Results.sites;
if ischar(bank)
    bank = str2double(bank);
end
if numel(bank) > 1
    bank = bank(:);
end
if size(implant,1) > 1
    error('You need to change code if you really put multiple Neuropixels probes in one rat')
end
is_neuropixels = contains(implant.implant_type, 'neuropixels', 'IgnoreCase', true) & ...
                  ~isnan(implant.neuropixels_sn);
if ~is_neuropixels
    error('This implant is not a Neuropixels probe.')
end
is3A=contains(implant.implant_type, '3A');
is3B=contains(implant.implant_type, '3B');
if (~is3A && ~is3B) || (is3A && is3B)
    error('Is this a 3A or 3B probe?')
end
if is3A
    % based on page 9 of the Neuropix Phase 3A End-user manual
    % The only probe that worked was option 1. Assume this is option 1
    length_without_contact_um = 137;
elseif is3B
    % page 9 of the Neuropixels phase 3B manual    
    length_without_contact_um = 195;
end
    
dist_um = length_without_contact_um + 10 + repmat(0:1:191,2,1) * 20;
dist_um = dist_um(:);

dist_um = dist_um(sites); 

dist_um = dist_um + bank * 20 * 192;