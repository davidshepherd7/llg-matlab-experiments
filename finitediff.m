function [y_approx, mesh_points] = finitediff(input_function,N,y0,y1)
% Solve -y''(x) = f(x) on (0,1) using finite difference method with N mesh
% points (including boundaries). Boundary conditions are set by y0 and y1. 
% input_function should be a function handle for f, might not work for
% constant f because of vectorisation problems.

h = 1/(N-1);
mesh_points = 0:h:1;

% hack to force vectorisation of constant input functions:
input_function2 = @(x) (input_function(x) + 0.*x);

% calculate f:
f = input_function2(mesh_points)';
f(1) = y0; f(end) = y1; % add boundary conditions to f vector

% construct matrix with 2 on diagonals, -1 on either side of diagonal then
% multiply by 1/h^2
r = [2,-1,zeros(1,N-2)];
A = (1/(h^2)) * toeplitz(r);

% Change top and bottom rows to contain a single 1, corresponding to the
% equation y(x_0) = y0, y(x_N) = y1
A(1,:) = zeros(1,N); A(N,:) = A(1);
A(1,1) = 1; A(N,N) = 1;

y_approx = A\f;

end