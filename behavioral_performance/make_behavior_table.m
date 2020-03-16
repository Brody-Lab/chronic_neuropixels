% MAKE_BEHAVIOR_TABLE Make a CSV table with information on the sessions
% used for measuring performance, including session ID, date, rat name,
% whether the rat was tethered for recording
function [] = make_behavior_table()
    add_folders_to_path;
    P = get_parameters;
    Rat_info = readtable(P.rat_info_path);
    T = struct;
    T.sessid=[];
    T.rat="";
    T.date = datetime.empty;
    Rec_sess = readtable(P.recording_sessions_path);
    for i = 1:numel(Rat_info.rat_name)
        dates_tethered = Rec_sess.date(strcmp(Rec_sess.rat, Rat_info.rat_name{i}));
        % add the sessids for the dates when the rat was tethered
        date_char = concatenate_for_sql(datetime(dates_tethered));
        [sessid, dates_tethered]=bdata(['select sessid, sessiondate from sessions where ratname="{S}" and sessiondate in (' date_char ...
                      ') order by sessiondate'], Rat_info.rat_name{i});
        T = add_to_table(T, dates_tethered, Rat_info.rat_name{i}, sessid, true);
        % add the sessids for the dates when the rat was in the training room
        date_surgery = datestr(datetime(Rat_info.Neuropixels_surgery_date(i)), 'yyyy-mm-dd');
        date_30d_b4 = datestr(datetime(Rat_info.Neuropixels_surgery_date(i))-30, 'yyyy-mm-dd');
        [dates_training, sessid] = bdata(['select sessiondate, sessid from sessions where sessiondate >="' date_30d_b4 ...
                                '" and sessiondate <= "' date_surgery '" and ratname="' Rat_info.rat_name{i} '"']);
        T = add_to_table(T, dates_training, Rat_info.rat_name{i}, sessid, false);
    end
    T = struct2table(T);
    writetable(T, P.behavior_table_path);
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
%   sessid
%       A numeric array specifying session ID's
%
%=OUTPUT
%   T
%       The tabular structure after new elements have been added
function T = add_to_table(T, dates, rat, sessid, tethered)
    dates = datetime(dates, 'format', 'yyyy-MM-dd');
    if numel(dates)~=numel(sessid)
        error('DATES and SESSID have different numbers of elements')
    end
    inds = numel(T.date)+1:numel(T.date)+numel(dates);
    T.rat(inds,1) = string(rat);
    T.date(inds,1) = dates;
    T.sessid(inds,1) = sessid;
    T.tethered(inds,1) = tethered;
end