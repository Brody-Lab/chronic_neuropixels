% GET_DATE_OF_IMPLANT return the date of an implant of one of Thomas's rat
%
%=INPUT
%
%   rat_name
%       A char, string, or cell array specifiying the name of a rat, e.g.
%       'T219'
%
%   probe_sn
%       The serial number of a probe
%
%=OUTPUT
%   
%   implant_date
%       A datetime element
function [implant_date] = get_date_of_implant(rat_name, probe_sn)
P =get_parameters;
Log = readtable(P.implant_log_path);
idx = strcmp(Log.rat, rat_name) & Log.neuropixels_sn == probe_sn;
if sum(idx) < 1
    error('Cannot find an entry in the implant log matching the %s and %i', rat_name, probe_sn)
elseif sum(idx) > 1
    error('Found multiple entries in the implant log matching the %s and %i', rat_name, probe_sn)
else
    implant_date = Log.implant_date(idx);
end