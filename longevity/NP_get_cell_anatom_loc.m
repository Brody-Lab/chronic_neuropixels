% NP_GET_CELL_ANATOM_LOC get anatomical positions of cells
%
%   Refer to the two following pages for specifying angles:
%       https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
%       https://www.mathworks.com/help/phased/ref/roty.html
%
%=INPUT
%
%   implant
%       a row of the imported Implants_log
%
%   dist_from_tip_um
%       A vector specifying the distance of each cell from the tip
%
%=OUTPUT
%   anatom_loc_mm
%       A structure with the AP, ML, and DV positions. A negative ML
%       value indicates that it is in the left hemisphere. A negative DV
%       value indicates that it is below the reference landmark 
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   reference
%       A char, string, or cell array specifying the reference landmark:
%       either 'bregma' or 'IA0';

function anatom_loc_mm = NP_get_cell_anatom_loc(implant, dist_from_tip_um, varargin)
validateattributes(implant, {'table', 'struct'},{'nonempty'})
parseobj = inputParser;
addParameter(parseobj, 'reference', 'bregma', @(x) ismember(x, {'bregma', 'IA0'}))
parse(parseobj, varargin{:});
P = parseobj.Results;
if implant.reference == "IA0" && P.reference == "bregma"
    implant.AP_mm = implant.AP_mm - 9;
    implant.DV_mm = implant.DV_mm - 10;
elseif implant.reference == "bregma" && P.reference == "IA0"
    implant.AP_mm = implant.AP_mm + 9;
    implant.DV_mm = implant.DV_mm + 10;
end
yaw = implant.horizontal_deg; % most of time is zero
pitch = implant.sagittal_deg;
roll = implant.coronal_deg;
R = rotz(yaw)*roty(pitch)*rotx(roll); % non-communicative
% made negative to be consistent with the right hand rule
u(3,:) = dist_from_tip_um(:)'/1000 - implant.depth_mm;
v = R*u + [implant.AP_mm; implant.ML_mm; implant.DV_mm];
anatom_loc_mm.AP = v(1,:)';
anatom_loc_mm.ML = v(2,:)';
anatom_loc_mm.DV = v(3,:)';