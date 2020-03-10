% MAKE_RECORDING_SESSIONS_TABLE Make a CSV specifying the recording
% sessions. Because this requires the "tzluo" repository, it is not called
% by any other function in the "chronic_neuropixels" repository
function [] = make_recording_sessions_table()
    P = get_parameters;
    if ~isfolder('Z:\RATTER\PhysData\Raw\')
        error('Please mount the ARCHIVE server')
    end
%     if ~isfolder(P.tzluo_path)
%         error('Please clone the TZLUO repository from the Brody Lab gitlub')
%     end
%     fprintf('Temporarily adding the TZLUO repository to the search path.\n')
%     addpath(genpath(P.tzluo_path))   
    T = struct;
    T.rat = "";
    T.date = datetime.empty;
    %K265
    T = add_to_table(T, 'K265', '2019-05-27');
    %T173
    T = add_to_table(T, 'T173', '2018-04-16');
    %T176
    Pharma = PB_import_pharmacology_log;
    idx = Pharma.Rat == 'T176' & Pharma.Manipulation == 'saline';
    T = add_to_table(T, 'T176', Pharma.Date(idx));
    %T181
    OE = PB_import_opto_ephys_log('T181');
    dates = OE.T181.Date;
    T = add_to_table(T, 'T181', dates);
    %T212
    dates = {'8/9/2019';...
             '8/13/2019';...
             '8/14/2019';...
             '8/15/2019';...
             '8/16/2019';...
             '8/18/2019';...
             '8/19/2019';...
             '8/20/2019';...
             '8/21/2019'};
     T = add_to_table(T, 'T212', dates);
     % T219
     T = add_opto_ephys_rat(T, 'T219');
     % T223
     T = add_opto_ephys_rat(T, 'T223');
     % T224
     dates = {'2020-01-03'};
     T = add_to_table(T, 'T224', dates);
     % T249
     dates = {'2020-02-11';...
             '2020-02-12'};
     T = add_to_table(T, 'T249', dates);
     % save
     writetable(struct2table(T), P.recording_sessions_path)
%      rmpath(genpath(P.tzluo_path))
%      fprintf('Removed the TZLUO repository to the search path.\n')
end
%% ADD_TO_TABLE
% add dates to a tabular structure
%=INPUT
%   T
%       A structure in a tabular form
%
%   rat_name
%       A char or string, does not have to match the number of elements of
%       DATES
%
%   dates
%       A DATETIME array
%
%=OUTPUT
%   T
%       The tabular structure after new elements have been added
function T = add_to_table(T, rat_name, dates)
    dates = datetime(dates, 'format', 'yyyy-MM-dd');
    k = numel(T.date);
    inds = k+1:k+numel(dates);
    T.rat(inds,1) = rat_name;
    T.date(inds,1) = dates;
end
%% ADD_OPTO_EPHYS_RAT
% add dates of a recording session of an opto-ephys rat
%=INPUT
%   T
%       A structure in a tabular form
%
%   rat_name
%       A char or string, does not have to match the number of elements of
%       DATES
%
%=OUTPUT
%   T
%       The tabular structure after new elements have been added
function T = add_opto_ephys_rat(T, rat_name)
    file_list = dir(['Z:\RATTER\PhysData\Raw\Thomas\' rat_name]);
    file_list=file_list(arrayfun(@(x) numel(x.name)>=15, file_list));
    dates = arrayfun(@(x) x.name(6:15), file_list, 'uni', 0);
    dates = datetime(dates, 'inputformat', 'yyyy_MM_dd');
    OE = PB_import_opto_ephys_log(rat_name);
    dates = intersect(dates, OE.(rat_name).Date);
    T = add_to_table(T, rat_name, dates);
end