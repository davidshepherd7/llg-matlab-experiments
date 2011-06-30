function [y_approx, mesh_points] = finitediff(input_function,N,y0,y1)
% Solve -y''(x) = f(x) on (0,1) using finite difference method with N mesh
% points (including boundaries). Boundary conditions are set by y0 and y1. 
% input_function should be a function handle for f, might not work for
% constant f because of vectorisation problems.

h = 1/(N-1)
mesh_points = 0:h:1

% hack to force vectorisation of constant input functions:
input_function = @(x) (input_function(x) + 0.*x);

% calculate f:
f = input_function(h:h:1-h)';
f(1) = f(1) + y0/(h^2);
f(end) = f(end) + y1/(h^2)

% construct matrix with 2 on diagonals, -1 on either side of diagonal then
% multiply by 1/h^2
r = [2,-1,zeros(1,N-4)];
A = (1/(h^2)) * toeplitz(r);

y_approx = A\f;

y_approx = [y0; y_approx; y1]

end