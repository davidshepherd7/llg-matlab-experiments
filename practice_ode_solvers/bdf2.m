function [Y, t] = bdf2(fcn, y0, h, N_steps)
% Solve the ODE Y'(t) = fcn(t,Y) using the BDF2 method. Initial value is
% calculated by ??. y_n+1 is calculated from the implicit formula at each
% step by fixed point iteration with the inital guess given by forward
% Euler.

% Issues:
% Stability seems to be limited by the fixed point iterations.
% RK method takes input function which is the transpose of my set up so
% needs a (wasteful) transpose to fix it


% Initialise
t = 0:h:h*(N_steps-1);
N_eqns = length(y0);
Y = zeros(N_steps,N_eqns);
implicit_y_accuracy = 0.005;

Y(1,:) = y0;
% Calculate second "initial value" using RK4 (built in - cheating a bit)
[unused Y2] = ode45(fcn, [0 h], y0);
Y(2,:) = Y2(end,:);

% Time stepping:
for i= 3:N_steps
    
    % Have y_n+1 implicitly so we need to solve for it:
    Y_temp = Y(i-1,:) + h*feval(fcn, t(i), Y(i-1,:))'; % use forward euler to get initial guess
    diff = 10*implicit_y_accuracy;   % ensure while loop is never accurate enough on first try
    j = 0;  % avoid infinite loops by counting number of iterations
    
    while (diff > implicit_y_accuracy)&&(j<10)
        Y_temp_prev = Y_temp;
        Y_temp = (2*h/3)*feval(fcn, t(i), Y_temp) + (4/3)*Y(i-1,:) - (1/3)*Y(i-2,:);    %BDF2 formula
        diff = max( abs(Y_temp(:) - Y_temp_prev(:)) );  % Measurement of change between steps
        j = j+1;
    end
    
    % Give a warning if Y did not converge
    if j>=10; warning('backward_euler:NCY','A value of Y did not converge in 10 or less steps'); end
    
    % Put into solution vector
    Y(i,:) = Y_temp;
    
end