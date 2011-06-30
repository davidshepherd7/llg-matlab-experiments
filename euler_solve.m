function [y_sol,t] = euler_solve(fcn, y0, h, N_steps)

% Solve ODE Y' = A*Y + fcn where A is a matrix and Y, fcn are vectors of
% functions, Y(t=0) = y0 using first order forward Euler method.
% Step size is given by h, number of steps is given by N_steps.

% Initialise vectors
N_eqns = length(y0);
y_sol = zeros(N_steps,N_eqns);
y_sol(1,:) = y0;

t = 0:h:h*(N_steps-1);


% Solve:
for i= 2:N_steps
    y_sol(i,:) = y_sol(i-1,:) + h*feval(fcn,t(i-1),y_sol(i-1,1),y_sol(i-1,2));
end