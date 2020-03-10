% IMPORT_OPTO_EPHYS_LOG A table of optogenetics experiments
%
%   LOG = IMPORT_OPTO_EPHYS_LOG() makes a structure with all data in Opto-ephys log.xlsx
%   LOG = IMPORT_OPTO_EPHYS_LOG(Rats) makes a table with only the rats specified in the string
%   or cell array "Rats"
%   LOG = IMPORT_OPTO_EPHYS_LOG(Rats, Dates) makes a table 
function [Log] = import_opto_ephys_log(varargin)
input_parser = inputParser;
addOptional(input_parser, 'rat_subset', '', @(x) ischar(x) || isstring(x) || iscell(x))
parse(input_parser, varargin{:});
for param = fields(input_parser.Results)'; eval([param{:} ' = input_parser.Results.' param{:} ';']); end
%% Import the spreadsheet
kFibers = {'AL', 'AR', 'PL', 'PR', 'sham'};
P = get_parameters;
warning('off', 'MATLAB:table:ModifiedVarnames');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
Log = struct;
[~, Rats] = xlsfinfo(P.opto_ephys_log_path); % the name of each sheet is the name of a rat;
Rats = string(Rats);
if ~isempty(rat_subset)
    if ischar(rat_subset)
        Rats = Rats(Rats == rat_subset);
    else
        Rats = Rats(any(Rats(:)' == rat_subset(:),1));
    end
end
for r = 1:numel(Rats)
    rat = Rats{r};
    Log.(rat) = readtable(P.opto_ephys_log_path, 'Sheet', rat);
    [~, Log.(rat).Date] = standardize_datestr(datetime(Log.(rat).Date));
end
warning('on', 'MATLAB:table:ModifiedVarnames');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
%% Delete rows without a date
for rat = Rats; rat = rat{:};
    Log.(rat)(~arrayfun(@isdatetime, Log.(rat).Date), :) = [];
    if isempty(Log.(rat))
        Log = rmfield(Log, rat);
    end
end
Rats = string(fields(Log))';
%% Reorder the dates
for rat = Rats; rat = rat{:};
    [~, Idx] = sort(Log.(rat).Date);
    Log.(rat) = Log.(rat)(Idx,:);
end
%% Change some columns from 1 & NaN to true/false
for rat = Rats; rat = rat{:};
    if isvar(Log.(rat), 'isAnesthetized')
        Log.(rat).isAnesthetized = ~isnan(Log.(rat).isAnesthetized);
    else
        Log.(rat).isAnesthetized = false(numel(Log.(rat).Date), 1);
    end
    if isvar(Log.(rat), 'noEphys')
        Log.(rat).noEphys = ~isnan(Log.(rat).noEphys);
    else
        Log.(rat).noEphys = true(numel(Log.(rat).Date), 1);
    end
    if isvar(Log.(rat), 'isConcerning') 
        if iscell(Log.(rat).isConcerning)
            Log.(rat).isConcerning = ~cellfun(@isempty, Log.(rat).isConcerning);
        else
            Log.(rat).isConcerning = ~isnan(Log.(rat).isConcerning);
        end
    else
        Log.(rat).isConcerning = false(size(Log.(rat).Date));
    end
end
%% Delete rows with isConcerning == true
for rat = Rats; rat = rat{:};
    Log.(rat) = Log.(rat)(~Log.(rat).isConcerning, :);
end
%% Add sess_id
% Rows with the same date are entered such that the later sessions are in lower rows
fprintf('\nChecking Opto-ephys log.xlsx against the MySQL table ...')
for rat = Rats; rat = rat{:};
    if ~isvar(Log.(rat), 'sessid')
        Log.(rat).sessid = nan(numel(Log.(rat).Date),1);
    end
    idx_to_remove = [];
    for d = 1:numel(Log.(rat).Date)
        if isvar(Log.(rat), 'misnamed_as') && ~isempty(Log.(rat).misnamed_as{d})
            ratname = Log.(rat).misnamed_as{d};
        else
            ratname = rat;
        end
        
        date_str = datestr(Log.(rat).Date(d), 'yyyy-mm-dd');
        command = sprintf(['select sessid, hostname from sessions where ratname = "%s" and sessiondate = "%s" ' ...
                          'order by starttime'], ratname, date_str);
        [sessid, hostname] = bdata(command);
        hostname = cellfun(@(x) str2double(x(4:end)), hostname);
        sessid = sessid(hostname > 200 | hostname == 19); % phys rigs only rigs 
        if numel(sessid) > 1
            if isnan(Log.(rat).sessid(d))
                Log.(rat).sessid(d) = resolve_mult_sess_per_day(ratname, Log.(rat).Date(d), 'Log', Log);
            else
                is_match = sessid == Log.(rat).sessid(d);
                if sum(is_match) ~= 1
                    Log.(rat).sessid(d) = resolve_mult_sess_per_day(ratname, Log.(rat).Date(d));
                end                
            end
        elseif numel(sessid) == 0
            warning(sprintf('%s: Could not find sessid for %s %s. Ignoring this row in the Opto-ephys log.', ...
                     mfilename, rat, Log.(rat).Date(d)))
            idx_to_remove = [idx_to_remove; d];
        else
            if isnan(Log.(rat).sessid(d))
                Log.(rat).sessid(d) = sessid;
            elseif sessid ~= Log.(rat).sessid(d)
                error('%s: For %s %s, the sessid in the Opto-ephys log is %i, but it is %i in the MySQL table', ...
                      mfilename, rat, Log.(rat).Date(d), Log.(rat).sessid(d), sessid)
            end
        end
    end
    Log.(rat)(idx_to_remove,:) = [];
end
fprintf(' done\n')
%% IsIllum - Which fiber was illuminated?
% Assume the set of possible fiber positions are {'AL', 'AR', 'PL', 'PR} for now, but this will likely change

% get the default encoding in number of bytes for unicode characters for each fiber
fiber_bytes = struct;
for fibers = kFibers; fib = fibers{:};
    fiber_bytes.(fib) = unicode2native(fib);
end

% Preset as false
IsIllum = struct;
for rat = Rats; rat = rat{:};
for fibers = kFibers; fib = fibers{:};
    nSess = numel(Log.(rat).Date);
    IsIllum.(rat).(fib) = false(nSess,1);    
end
end

% parse the text
for rats = Rats; rat = rats{:};
for s = 1:numel(Log.(rat).Date) % number of sessions
for fibers = kFibers; fib = fibers{:};        
    if isvar(Log.(rat), 'FiberIlluminated') && ...
       ~isempty(regexp(Log.(rat).FiberIlluminated{s}, fib, 'once'))
        IsIllum.(rat).(fib)(s,1) = true;            
    end
end
end
end

% *********
% * Check *
% *********
% Are there sessions with a non-useful combination of illuminations?
for rats = Rats; rat = rats{:};
    idx_confounded = IsIllum.(rat).AL & IsIllum.(rat).PR | ...
                     IsIllum.(rat).PL & IsIllum.(rat).AR;
    if any(idx_confounded)
        error('You stimulated different brain areas in different hemispheres. Sir, you better have a good reason.')
    end
end

% determine whether the inactivation was unilateral or bilateral 
%   0 = left hemisphere
%   1 = right hemisphere
%   2 = bilateral
% It takes care of sessions in which both cortex and striatum were stimulated
for rats = Rats; rat = rats{:};
    if isvar(Log.(rat), 'FiberIlluminated')
        n_sess = numel(Log.(rat).FiberIlluminated);
        IsIllum.(rat).laterality = strings([n_sess,1]);

        idx_bilat =  IsIllum.(rat).AL &  IsIllum.(rat).AR | ...
                     IsIllum.(rat).PL &  IsIllum.(rat).PR;
        idx_left  =  IsIllum.(rat).AL & ~IsIllum.(rat).AR | ...
                     IsIllum.(rat).PL & ~IsIllum.(rat).PR;
        idx_right = ~IsIllum.(rat).AL &  IsIllum.(rat).AR | ...
                    ~IsIllum.(rat).PL &  IsIllum.(rat).PR;

        IsIllum.(rat).laterality(idx_left, 1) = string('L');
        IsIllum.(rat).laterality(idx_right,1) = string('R');
        IsIllum.(rat).laterality(idx_bilat,1) = string('B');
    end
end
%% Prune rats
if nargin > 0 && ~isempty(varargin{1})
    Rats_to_include = string(varargin{1}); % strings are nice to work with

    % *********
    % * Check *
    % *********
    % Is each rat that the user wants included within the spreadsheet    
    for r = 1:numel(Rats_to_include)
        if ~any(Rats_to_include{r} == Rats)
            warning([Rats_to_include{r} ' is not a rat in the Opto-ephys log'])
        end
    end
    
    % Remove rats that aren't on the list to be included
    for rats = Rats; rat = rats{:};
        if ~any(rat == Rats_to_include)
            Log = rmfield(Log, rat);
        end
    end
end
%% Prune dates
if nargin > 1 && ~isempty(varargin{2})
    Dates_to_include = datetime(varargin{2}); % datetimes are nice to work with
    
    for rats = Rats; rat = rats{:};
        if ~isvar(Log, rat); continue; end % might have been removed in the previous section
        [idx] = ismember(Log.(rat).Date, Dates_to_include);
        Log.(rat) = Log.(rat)(idx,:);
    end
end
%% right_power_multiplier
for rats = Rats; rat = rats{:};
    if isvar(Log.(rat), 'right_power_multiplier')
        Log.(rat).right_power_multiplier(isnan(Log.(rat).right_power_multiplier)) = 1;
    else
        Log.(rat).right_power_multiplier = nan(numel(Log.(rat).Date), 1);
    end
end
%% Remove "isConcerning" sessions
for rats = Rats; rat = rats{:};
    Log.(rat) = Log.(rat)(~Log.(rat).isConcerning, :);
end