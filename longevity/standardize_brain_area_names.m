function Cells = standardize_brain_area_names(Cells)
% STANDARDIZE_BRAIN_AREA_NAMES Assign standardized names to each unit and
% each electrode
%
%=INPUT
%
%   Cells
%       A structure made using "collect_cells_files" and
%       "postprocess_cells"
%
%=OUTPUT
%
%   Cells
%       The input structure with the fields "region_names" and
%       "electrodes.brain_area" modified
%

for i = 1:numel(Cells)
    rat_name = unique(string(Cells{i}.rat));
    if numel(rat_name) > 1
        error('Multiple unique rat names in a Cells file.')
    end
    switch rat_name
        case {'T173', 'T181', 'T182'}
            Cells{i}.region_names(Cells{i}.DV >  -2.4) = "MCtx";
            Cells{i}.region_names(Cells{i}.DV <= -2.4 & ...
                                  Cells{i}.DV >  -5) = "dmStr";
            Cells{i}.region_names(Cells{i}.DV <= -5 & ...
                                  Cells{i}.ML <= 3.1) = "NAc";
            Cells{i}.region_names(Cells{i}.ML > 3.1) = "piriform";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV >  -2.4) = "MCtx";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV <= -2.4 & ...
                                           Cells{i}.electrodes.DV >  -5) = "dmStr";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV <= -5 & ...
                                           Cells{i}.electrodes.ML <= 3.1) = "NAc";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.ML >  3.1) = "piriform";
        case {'T176', 'T212', 'T224', 'T249'}
            Cells{i}.region_names(Cells{i}.region_names=="M2" | ...
                                  Cells{i}.region_names=="Cg1") = "dmFC";
            Cells{i}.region_names(Cells{i}.region_names=="PrL" | ...
                                  Cells{i}.region_names=="MO") = "vmFC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="M2" | ...
                                           Cells{i}.electrodes.brain_area=="Cg1") = "dmFC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="PrL" | ...
                                           Cells{i}.electrodes.brain_area=="MO") = "vmFC";
        case {'T219', 'T223', 'T227'}
            Cells{i}.region_names(Cells{i}.region_names=="M1") = "MCtx";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="M1") = "MCtx";
            Cells{i}.region_names(Cells{i}.region_names=="dStr") = "dmStr";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="dStr") = "dmStr";
            Cells{i}.region_names(Cells{i}.region_names=="vStr") = "NAc";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="vStr") = "NAc";
        case 'T170'
            Cells{i}.region_names(Cells{i}.DV >   -1.6) = "dmFC";
            Cells{i}.region_names(Cells{i}.DV <=  -1.6) = "vmFC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV >   -1.6) = "dmFC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV <=  -1.6) = "vmFC";
        case 'K265'    
            Cells{i}.region_names(Cells{i}.DV >  -2.4) = "MCtx";
            Cells{i}.region_names(Cells{i}.DV <= -2.4 & ...
                                  Cells{i}.DV >  -5) = "dmStr";
            Cells{i}.region_names(Cells{i}.DV <= -5) = "NAc";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV >  -2.4) = "MCtx";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV <= -2.4 & ...
                                  Cells{i}.electrodes.DV >  -5) = "dmStr";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV <= -5) = "NAc";            
        case 'T179'
            Cells{i}.region_names(Cells{i}.region_names=="ECIC" | ...
                                  Cells{i}.region_names=="MiTG") = "nIC";
            Cells{i}.region_names(Cells{i}.region_names=="Post") = "postsubiculum";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="ECIC" | ...
                                           Cells{i}.electrodes.brain_area=="MiTG") = "nIC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="Post")="postsubiculum";
            Cells{i}.region_names(Cells{i}.region_names=="V1M" | ...
                                  Cells{i}.region_names=="RSGa") = "dpCtx";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="V1M" | ...
                                           Cells{i}.electrodes.brain_area=="RSGa") = "dpCtx";
        case {'T196', 'T209'}
            Cells{i}.region_names(Cells{i}.region_names=="IC" | ...
                                  Cells{i}.region_names=="aIC") = "nIC";
            Cells{i}.region_names(Cells{i}.region_names=="sSC" | ...
                                  Cells{i}.region_names=="iSC" | ...
                                  Cells{i}.region_names=="dSC") = "SC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="IC" | ...
                                           Cells{i}.electrodes.brain_area=="aIC") = "nIC";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="sSC" | ...
                                           Cells{i}.electrodes.brain_area=="iSC" | ...
                                           Cells{i}.electrodes.brain_area=="dSC") = "SC";
        case 'A230'
            Cells{i}.region_names(Cells{i}.region_names=="DMS") = "dmStr";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="DMS") = "dmStr";            
        case 'A241'
            if Cells{i}.probe_serial=="18194823302"
                Cells{i}.region_names(Cells{i}.region_names=="DMS") = "dmStr";   
                Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="DMS") = "dmStr";                            
            elseif Cells{i}.probe_serial=="18194823631"
                Cells{i}.region_names(Cells{i}.region_names=="dorsomedial striatum") = "dmStr";   
                Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="dorsomedial striatum") = "dmStr";                            
                Cells{i}.region_names(Cells{i}.region_names=="bed nucleus of the stria terminalis") = "BNST";     
                Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="bed nucleus of the stria terminalis") = "BNST";                            
                Cells{i}.region_names(Cells{i}.region_names=="M1") = "MCtx";   
                Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="M1") = "MCtx";                            
            else
                error('Unknown probe serial number')
            end
        case 'A242'
            Cells{i}.region_names(Cells{i}.region_names=="lateral amygdala") = "amygdala";
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="lateral amygdala") = "amygdala";
        case 'A243'
            Cells{i}.region_names(Cells{i}.region_names=="dorsomedial striatum") = "dmStr";    
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="dorsomedial striatum") = "dmStr";       
            Cells{i}.region_names(Cells{i}.region_names=="M1") = "MCtx";   
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="M1") = "MCtx";                
        case 'A249'
            Cells{i}.region_names(Cells{i}.region_names=="ADS") = "dmStr";  
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="ADS") = "dmStr";                                        
            Cells{i}.region_names(Cells{i}.DV<-5.3) = "NAc";     
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.DV<-5.3) = "NAc";                           
            Cells{i}.region_names(Cells{i}.region_names=="M2") = "MCtx";            
            Cells{i}.electrodes.brain_area(Cells{i}.electrodes.brain_area=="M2") = "MCtx";                                        
        otherwise
            error('Unknown rat')
    end
end