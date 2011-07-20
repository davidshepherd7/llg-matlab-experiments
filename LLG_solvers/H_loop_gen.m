function [H_loop] = H_loop_gen(H_min,H_max,N,theta)
% Generate a loop of fields (in the z direction) from H_min to H_max and
% back again in N_steps steps. The field is at angle theta to the z
% axis.

% Set default theta =0
if nargin < 4; theta = 0; end

% Generate loop in field strength
H_loop_r1 = [ H_min * ones(N/(2*10),1) ; zeros(4*N/(2*10),1); H_max * ones(5*N/(2*10),1)]; % min to max
H_loop_r2 = flipud(H_loop_r1);                                                 % and back

H_loop_r = [ H_loop_r1; H_loop_r2];

% Align applied field at angle theta to z axis (phi = 0)
% [Using SI sph polar notation (r, theta=polar/inclination, phi=azimuthal)]
% [Note that output is in cartesian coordinates]
H_loop = zeros(N,3);
for i=1:N
    H_loop(i,:) = sphtocart([H_loop_r(i),theta,0]);
end