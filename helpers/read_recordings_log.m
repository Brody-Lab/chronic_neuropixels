function recordings_table = read_recordings_log(csvfile)
    %% imports a single sheet from the Neuropixels recordings_log Google sheet (downloaded as a local .csv file)
    %% outputs a Matlab table
    opts = detectImportOptions(csvfile,'TextType','string');
    opts = opts.setvartype('date','datetime');
    opts=opts.setvaropts('date','InputFormat','yyyy_MM_dd');    
    opts = opts.setvaropts('date','DatetimeFormat','dd-MMM-uuuu');    
    recordings_table=readtable(csvfile,opts);
end