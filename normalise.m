function [y_out] = normalise(y,Ms)

if nargin <2; Ms = 1; end

y_out = Ms*(y/norm(y,2));