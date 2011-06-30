function Y_rotated = rotate(Y,theta,phi)
% Rotate vector Y by theta and phi (Where theta is the inclination from z,
% phi is the angle from x. Both angles are measured clockwise.)

% Rotation matrices from Wolfram
R_phi = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];

R_theta = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];

Y_rotated = R_phi * (R_theta * Y);