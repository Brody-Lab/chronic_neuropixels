function age = get_age_at_implant(rat_name, probe_serial)
% GET_AGE_AT_IMPLANT return the age of the animal at days for an implant
% specified by the name of the animal and the serial number of the probe
%
%=INPUT
%
%   rat_name
%       A char vector, string, or cell string of the name of a single rat
%
%   probe_serial
%       A numeric specifying the serial number
%
%=OUTPUT
%
%   age
%       age of the animal on the day of the implant (in days from when the
%       animal was delivered to the animal facility, not from birth).

P = get_parameters;
T = readtable(P.age_at_implant_path);

idx = T.rat == rat_name & T.probe_serial == probe_serial;
if sum(idx)~=1
    error('Cannot find a unique age for animal %s and probe serial %i', rat_name, probe_serial)
end
age = T.age_day(idx);