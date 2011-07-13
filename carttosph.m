function X_sph =  carttosph(X_cart)
% Using ISO standard coordinate system:
% [r, θ, φ]	[x, y, z] = r [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
% Results in radians.
% Note that default MATLAB commands use a different system!

% r = sqrt(x^2 + y^2 + z^2)
X_sph(1) = sqrt(X_cart(1)^2 + X_cart(2)^2 + X_cart(3)^2);

% θ = acos(z/r)
X_sph(2) = acos(X_cart(3)/X_sph(1));

% φ = atan(y/x)
X_sph(3) = atan2(X_cart(2),X_cart(1));
