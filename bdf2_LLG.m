% Based on function:
%function [Y, t] = bdf2(LLG, y0, h, N_steps)
% Solve the ODE Y'(t) = LLG(t,Y) using the BDF2 method. Initial value is
% calculated by built in RK4 method. y_n+1 is calculated from the implicit formula at each
% step by fixed point iteration with the inital guess given by forward
% Euler.


% Notes: Stability seems to be limited by the fixed point iterations.
% Everything is normalised prior to storage in Y, this allows for
% new computations of H to be done from Y whenever they are needed.
% The function "normalise(y)" normalises y s.t. sqrt(y_1^2 +...) = 1.

%=========================================================================
% Input:
%=========================================================================
clear;

h = 0.1;
N_steps = 1000;
y0 = [1 0 0]';
N_eqns = length(y0);

% Magnetic parameters:
Ms = 1;
H_app = [0 0 0]';

% Ellipsoid geometry:
a = 2;  % long axis length
b = 1;  % short axis length
ellipsoid_theta = 0; %long axis elevation
ellipsoid_phi = 0;   %long axis azimuthal?

H_demag = ellipsoid_demag(a, b, y0,ellipsoid_theta,ellipsoid_phi);
H_demag_out = zeros(N_eqns,N_steps);
H_demag_out(:,1) = H_demag; % Store for use in output (hard to use properly since H is in the LLG function)
H = H_app + H_demag;

% LLG function:
gamma = 0.1;    % Gamma just affects the size of dM/dt. 
alpha = 0.7;    % The relative strength of damping is propotional to alpha.
LLG = @(t,M,H) (gamma/(1+alpha^2))*(cross(M,H)...
    - (alpha/Ms)*cross(M,cross(M,H)));
% built in methods use opposite type of vectors to mine so take transpose
LLG_transpose = @(t,M,H) (gamma/(1+alpha^2))*(cross(M,H)...
    - (alpha/Ms)*cross(M,cross(M,H)))';

%=========================================================================


% Initialise
t = 0:h:h*(N_steps-1);

Y = zeros(N_eqns,N_steps);
implicit_y_accuracy = 0.005;

Y(:,1) = normalise(y0);
% Calculate second "initial value" using a single Euler step
Y(:,2) = Y(:,1) + h*feval(LLG,t(1),Y(:,1),H);


% Main loop:
for i= 3:N_steps
    
    % In full solution solve for new H vector here
    H_demag = ellipsoid_demag(a, b, Y(:,i-1), ellipsoid_theta, ellipsoid_phi); % recalculate H_demag with new M
    H = H_app + H_demag;
    H_demag_out(:,i-1) = H_demag;
    
    % Have y_n+1 implicitly so we need to solve for it:
    Y_temp = normalise( Y(:,i-1) + h*feval(LLG, t(i), Y(:,i-1), H) ); % use forward euler to get initial guess
    diff = 10*implicit_y_accuracy;   % ensure while loop is run at least once
    j = 0;  % avoid infinite loops by counting number of iterations
    
    while (diff > implicit_y_accuracy)&&(j<10)
        Y_temp_prev = Y_temp;
        Y_temp = normalise( (2*h/3)*feval(LLG, t(i), Y_temp, H) + (4/3)*Y(:,i-1) - (1/3)*Y(:,i-2) );   %BDF2 formula
        diff = max( abs(Y_temp(:) - Y_temp_prev(:)) );  % Measurement of change between steps
        j = j+1;
    end
    
    % Give a warning if Y did not converge
    if j>=10; warning('backward_euler:NCY','A value of Y did not converge in 10 or less steps'); end
    
    % Put into solution vector
    Y(:,i) = Y_temp;
    
    % Use MATLAB's solver for comparison
    
end

%=========================================================================
% Output:
%=========================================================================

% plot line showing path of M (parameterised by t)
plot3(Y(1,:),Y(2,:),Y(3,:),0,0,0,'x')
axis([-1,1,-1,1,-1,1])
legend('Path of magnetisation with time', 'Origin','Location','NorthEast' )

% also plot each component of M against t
figure
plot(t,Y)
legend('Mx','My','Mz')

% also plot H_demagOut with time
figure
plot(t,H_demag_out)
legend('H demag x','H demag y','H demag z')