% MAKE_BEHAVIOR_TABLE 
function [] = make_behavior_table()
    add_folders_to_path;
    P = get_parameters;
    Rat_info = readtable(P.rat_info_path);
    T = struct;
    T.sessid=[];
    T.rat="";
    Rec_sess = readtable(P.recording_sessions_path);
    for i = 1:numel(Rat_info.rat_name)
        dates_tethered = [];
        switch Rat_info.rat_name{i}
            case 'A230'
                dates_tethered = bdata(['select sessiondate from sessions where '...
                                        'sessiondate>="2019-07-15" and '...
                                        'ratname="A230" and hostname="Rig205"']);
            case 'A242'
                sessid = bdata(['select sessid from sessions where '...
                                'sessiondate>="2019-06-18" and '...
                                'ratname="A242" and hostname="Rig205"']); 
                % for these sessions there is no physdata collected and video
                % shows rat not plugged in or video is broken so best to
                % exclude these sessions
                A242_in_physrig_but_maybe_not_tethered = [704440, ...
                                                          704117, ...
                                                          705282, ...
                                                          709379];
                sessid = setdiff(sessid, A242_in_physrig_but_maybe_not_tethered);
                dates_tethered = bdata(['select sessiondate from sessions s where s.sessid in (' ...
                                        concatenate_for_sql(sessid) ')']);
            otherwise
                dates_tethered = Rec_sess.date(strcmp(Rec_sess.rat, Rat_info.rat_name{i}));
        end
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

function T = add_to_table(T, dates, rat, sessid, tethered)
    dates = datetime(dates, 'format', 'yyyy-MM-dd');
    inds = numel(T.sessid)+1:numel(T.sessid)+numel(sessid);
    T.rat(inds,1) = string(rat);
    T.date(inds,1) = dates;
    T.sessid(inds,1) = sessid;
    T.tethered(inds,1) = tethered;
end