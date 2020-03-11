% GET_GAIN_NOISE_DATA % save the data from gain and noise measurements as
% CSV files
function [] = get_gain_noise_data(varargin)
    parseobj = inputParser;
    addOptional(parseobj, 'remake', false, @(x) islogical(x) && isscalar(x))
    parse(parseobj, varargin{:});
    P_input = parseobj.Results;
    P=get_parameters;
    Log = readtable(P.gain_noise_log_path);
    for i = 1:size(Log,1)
        recording_fldr_path = fullfile(P.gain_noise_data_path, Log.recording_id{i});
        csv_file_path = [P.gain_noise_fldr_path filesep Log.recording_id{i} '.csv'];
        if isfile(csv_file_path) && ~P_input.remake
            continue
        end
        Res = analyze_gain_noise(recording_fldr_path);
        % find the largest gain amplification for which measurements were made
        ind= find(all(~isnan(Res.ap.gain) & ~isnan(Res.ap.noise_uV)), 1, 'last');
        if isempty(ind)
            error('Cannot determine the input gain for which the data were recorded')
        end
        T =struct;
        T.gain = Res.ap.gain(:,ind);
        T.noise_uV = Res.ap.noise_uV(:,ind);
        T.amplification = repmat(Res.amplification(ind), numel(T.gain), 1);
        T = struct2table(T);
        writetable(T, csv_file_path)
    end
end