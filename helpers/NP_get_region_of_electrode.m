function [region_of_electrodes, Region_from_tip_um] = NP_get_region_of_electrode(implant, varargin)
input_parser = inputParser;
addParameter(input_parser, 'sites', 1:384, @(x) isnumeric(x))    % electrode index
addParameter(input_parser, 'bank', 0, @(x) isnumeric(x))
parse(input_parser, varargin{:});
sites = input_parser.Results.sites;
bank = input_parser.Results.bank;

rat = implant.rat;
if strcmp(rat, 'T182'); rat = 'T181'; end % T181 and T182 have their probes in the same parts of their brain

switch rat
    case {'T176', 'T212'}
        Region_from_tip_um.MO  = [0 393.92];
        Region_from_tip_um.PrL = max(Region_from_tip_um.MO)  + [0, 1891.6];
        Region_from_tip_um.Cg1 = max(Region_from_tip_um.PrL) + [0, 793.95];
        Region_from_tip_um.M2  = max(Region_from_tip_um.Cg1) + [0, 920.81];
    case {'T170', 'T173', 'T181', 'T182'}
        Region_from_tip_um.Pir = 137 + [0 208];
        Region_from_tip_um.En  = max(Region_from_tip_um.Pir)  + [0, 755];
        Region_from_tip_um.Str = max(Region_from_tip_um.En) + [0, 3718];
        Region_from_tip_um.cc  = max(Region_from_tip_um.Str) + [0, 876];
        Region_from_tip_um.M2  = max(Region_from_tip_um.cc) + [0, 2174];
    case 'T179'
        Region_from_tip_um.MiTG  = 137 + [0 443.8795];
        Region_from_tip_um.ECIC = max(Region_from_tip_um.MiTG)  + [0, 2214.4];
        Region_from_tip_um.Post = max(Region_from_tip_um.ECIC) + [0, 1122.8];
        Region_from_tip_um.RSGa = max(Region_from_tip_um.Post) + [0, 264.65];
        Region_from_tip_um.V1M  = max(Region_from_tip_um.RSGa) + [0, 1717.3];
    case 'T196'
        Region_from_tip_um.IC  = 137 + [0 1496.1];
        Region_from_tip_um.dSC = max(Region_from_tip_um.IC)  + [0, 1471.4];
        Region_from_tip_um.iSC = max(Region_from_tip_um.dSC) + [0, 471.4];
        Region_from_tip_um.sSC = max(Region_from_tip_um.iSC) + [0, 534.1];
        Region_from_tip_um.PiSt  = max(Region_from_tip_um.sSC) + [0, 845.5];
        Region_from_tip_um.RSG  = max(Region_from_tip_um.PiSt) + [0, 1170.2];
        Region_from_tip_um.V2  = max(Region_from_tip_um.RSG) + [0, 436.6];
    case 'K265'
        Region_from_tip_um.Str  = 137 + [0 3504.5];
        Region_from_tip_um.CC = max(Region_from_tip_um.Str)  + [0, 1001.2];
        Region_from_tip_um.M2 = max(Region_from_tip_um.CC) + [0, 2234.7];
    case 'T209'
        Region_from_tip_um.aIC = 137 + [0 2000];
        Region_from_tip_um.SC = Region_from_tip_um.aIC + [2000 3830];
    case 'T209'
        Region_from_tip_um.aIC = 137 + [0 2000];
        Region_from_tip_um.SC = Region_from_tip_um.aIC + [2000 3830];
    case {'T219', 'T223'}
        Region_from_tip_um.vStr = 137 + [0 2610];
        Region_from_tip_um.dStr = max(Region_from_tip_um.vStr) + [0 2710];
        Region_from_tip_um.M1 = max(Region_from_tip_um.dStr) + [0 2580];
    otherwise
        error('This rat has not been indicated to have a Neuropixels probe')
end
%% Region_of_elec
elec_dist_from_tip_um = NP_get_distance_from_tip(implant, 'sites', sites, 'bank', bank);

n_electrodes = numel(elec_dist_from_tip_um);

region_of_electrodes(n_electrodes,1)="";

the_brain_areas = fields(Region_from_tip_um)';
for brain_areas = the_brain_areas; area = brain_areas{:};
    idx = elec_dist_from_tip_um >= Region_from_tip_um.(area)(1) & ...
          elec_dist_from_tip_um <  Region_from_tip_um.(area)(2);
      
    region_of_electrodes(idx,1) = string(area);
end
