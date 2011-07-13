function [Y_soln, t] = backward_euler(fcn, y0, h, N_steps)
% Solve the differential equation Y' = fcn(t,Y) with initial value y0
% using the backwards Euler method. Step size = h, number of steps taken =
% N_steps. Implicit equation is solved by using forward Euler for an initial
% value in a fixed point iteration method.

% Other parameters
implicit_y_accuracy = 0.01; % accuracy to solve for y_n+1 to

% Initialise
N_eqns = length(y0);
Y_soln = zeros(N_steps,N_eqns);
Y_soln(1,:) = y0;
t = 0:h:h*(N_steps-1);

for i= 2:N_steps
    
    % Have y_n+1 implicitly so we need to solve for it:
    Y_temp = Y_soln(i-1,:) + h*feval(fcn, t(i), Y_soln(i-1,:)); % use forward euler to get initial guess
    diff = 10*implicit_y_accuracy;   % ensure while loop is never accurate enough on first try
    j = 0;  % avoid infinite loops by counting number of iterations
    
    while (diff > implicit_y_accuracy)&&(j<10)
        Y_temp_prev = Y_temp;
        Y_temp =  Y_soln(i-1,:) + h*feval(fcn, t(i), Y_temp);
        diff = max( abs(Y_temp(:) - Y_temp_prev(:)) );
        j = j+1;
    end
    
    % Give a warning if Y did not converge
    if j>=10; warning('backward_euler:NCY','A value of Y did not converge in 10 or less steps'); end
    
    % Put into solution vector
    Y_soln(i,:) = Y_temp;
end