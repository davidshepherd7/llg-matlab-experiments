function X_sph =  carttosph(X_cart,system)

% % Default to ISO system
% if (nargin <2)
%     system = 'ISO';
% end

% Using ISO standard coordinate system:
% r = radial, θ = inclination/polar, φ = azimuthal
% [r, θ, φ]	[x, y, z] = r [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
% Results in radians.
% Note that default MATLAB commands use a different system!

% r = sqrt(x^2 + y^2 + z^2)
X_sph(1) = sqrt(X_cart(1)^2 + X_cart(2)^2 + X_cart(3)^2);

% θ = acos(z/r)
X_sph(2) = acos(X_cart(3)/X_sph(1));

% φ = atan(y/x)
X_sph(3) = atan2(X_cart(2),X_cart(1));

% % Convert other systems to ISO:
% % System with the meanings of theta and phi switched:
% if strcmp(system,'swap_angles')
%     temp = X_sph(2);
%     X_sph(2) = X_sph(3);
%     X_sph(3) = temp;
%     
% % System with theta as the elevation:
% elseif (strcmp(system,'elevation'))
%     X_sph(2) = pi/2 - X_sph(2);
%     
% % Change azimuthal direction to be towards +ve y
% elseif (strcmp(system,'posi_azi'))
%     X_sph(3) = - X_sph(3);
% 
% elseif (~strcmp(system,'ISO'))
%     error('sphtocart:unknowncoordsystem','Chosen coordinate system does not exist, ensure argument is a string.');
% end