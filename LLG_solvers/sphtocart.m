function X_cart = sphtocart(X_sph,system)

% % Default to ISO system
% if (nargin <2)
%     system = 'ISO';
% end
% 
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
% % For unrecognised 'system' strings give an error
% elseif (~strcmp(system,'ISO'))
%     error('sphtocart:unknowncoordsystem','Chosen coordinate system does not exist, ensure argument is a string.');
% end


% Using ISO standard coordinate system:
% r = radial, θ = inclination/polar, φ = azimuthal
% [r, θ, φ]	[x, y, z] = r [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
% Note that default MATLAB commands use a different syste

% x = r sin(θ)cos(φ)
X_cart(1) = X_sph(1)*sin(X_sph(2))*cos(X_sph(3));

% y = r sin(θ)sin(φ)
X_cart(2) = X_sph(1)*sin(X_sph(2))*sin(X_sph(3));

% z = r cos(θ)
X_cart(3) = X_sph(1)*cos(X_sph(2));
