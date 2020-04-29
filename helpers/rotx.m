function rotmat = rotx(alpha)
%rotx     Rotation matrix around x-axis
%   ROTMAT = rotx(ALPHA) returns the rotation matrix, ROTMAT, that rotates
%   a point around the x-axis for an angle ALPHA (in degrees). The point is
%   specified in the form of [x;y;z], with the x, y, and z axes forming a
%   right-handed Cartesian coordinate system. With the x axis pointing
%   towards the observer, ALPHA is measured counter-clockwise in the y-z
%   plane.
%
%   ROTMAT is a 3x3 matrix. The rotation of the point can be achieved by
%   left-multiplying ROTMAT with the point's coordinate vector [x;y;z].
%
%   % Example:
%   %   Rotate a point, (0,1,0), around x-axis 45 degrees
%   %   counter-clockwise.
%
%   p = [0;1;0];
%   p = rotx(45)*p
    validateattributes(alpha,{'numeric'},{'nonnan','scalar','finite'},'rotx','alpha',1);
    % rotate in the direction of y->z, counter-clockwise
    rotmat = [1 0 0;0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];
end
