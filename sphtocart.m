function X_cart = sphtocart(X_sph)
% Using ISO standard coordinate system:
% [r, θ, φ]	[x, y, z] = r [sin(θ)cos(φ), sin(θ)sin(φ), cos(θ)]
% Note that default MATLAB commands use a different system!

% x = r sin(θ)cos(φ)
X_cart(1) = X_sph(1)*sin(X_sph(2))*cos(X_sph(3));

% y = r sin(θ)sin(φ)
X_cart(2) = X_sph(1)*sin(X_sph(2))*sin(X_sph(3));

% z = r cos(θ)
X_cart(3) = X_sph(1)*cos(X_sph(2));