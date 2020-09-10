function Cells = add_SNR_data(Cells)
    P = get_parameters;
    for i=65:length(Cells)
        fprintf('\n---------\nWorking on cells file %g of %g.\n---------\n',i,length(Cells));
        if isfield(Cells{i},'jrc_file')
            jrc_file = Cells{i}.jrc_file;
        else
            jrc_file = Cells{i}.meta.path.spikesorted_output;
        end
        if ~exist(jrc_file,'file')
            jrc_file = find_jrc_file(jrc_file);
        end
        if isfield(Cells{i},'rec')
            meta = Cells{i}.rec.ap_meta;
        else
            meta = Cells{i}.meta.ap_meta;
        end
        [parent,~,~] = fileparts(jrc_file);
        full_file = dir([parent,'\*full.prm']);
        if ~isempty(full_file)
            prm_path = fullfile(full_file(1).folder,full_file(1).name);
        else
            prm_file = dir([parent,'\*prm']);       
            if isempty(prm_file)
                warning('Could not locate a prm file. Assuming default qq factor of 5.');
                prm_path = '';
            else
                prm_path = fullfile(prm_file(1).folder,prm_file(1).name);
            end
        end
        fprintf('Loading %s ... ',jrc_file);tic;
        JRC = load(jrc_file);
        fprintf(' took %s.\n',timestr(toc));
        if numel(JRC.unitCount)~=numel(Cells{i}.ks_good)
            error('Mismatched unit numbers in JRC and Cells files.');
        end        
        if isempty(prm_path)
            Cells{i}.SNR = get_SNR_stats(JRC,meta,[],true);
        else
            Cells{i}.SNR = get_SNR_stats(JRC,meta,prm_path,true) ;       
        end
    end
end

function jrc_file_out = find_jrc_file(jrc_file)
    experimenters = {'Thomas','Adrian'};
    prefix = fullfile('X:','RATTER','PhysData','NP_sorted');
    if ~isdir(prefix)
        error('Default replacement prefix isn''t a path.');
    end
    for i=1:length(experimenters)
        where = strfind(jrc_file,experimenters{i});
        if ~isempty(where)
            jrc_file_out = fullfile(prefix,jrc_file(where:end));
            if exist(jrc_file_out,'file')
                fprintf('Found %s \n to replace %s.\n',jrc_file_out,jrc_file);
                return
            else
                warndlg(sprintf('Could not find attempted alternative file path:\n%s.\nYou''ll have to locate one manually.',jrc_file_out));
                [jrc_file_out,path] = uigetfile('*.mat','Locate jrc file manually',prefix);             
                jrc_file_out = fullfile(path,jrc_file_out);
            end
        end
    end

end

