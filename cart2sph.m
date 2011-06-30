function [vec_sph] = cart2sph(vec_cart)
%CART2SPH Transform Cartesian to spherical coordinates.
%   Overloaded version to take a single vector of co-ordinates

%   [TH,PHI,R] = CART2SPH(X,Y,Z) transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y,Z to spherical
%   coordinates (azimuth TH, elevation PHI, and radius R).  The arrays
%   X,Y, and Z must be the same size (or any of them can be scalar).
%   TH and PHI are returned in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
%   Class support for inputs X,Y,Z:
%      float: double, single
%
%   See also CART2POL, SPH2CART, POL2CART.

%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.12.4.2 $  $Date: 2005/05/31 16:31:21 $

x = vec_cart(1);
y = vec_cart(2);
z = vec_cart(3);

hypotxy = hypot(x,y);
r = hypot(hypotxy,z);
elev = atan2(z,hypotxy);
az = atan2(y,x);

vec_sph = [az, elev, r];