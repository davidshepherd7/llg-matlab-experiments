function [T_out, Y_out] = bdf2(odefun,tspan, y0, h, jacobian)
% Solve the ODE Y'(t) = odefun(t,Y) using the BDF2 method. Initial value is
% calculated by ??. y_n+1 is calculated from the implicit formula at each
% step by fixed point iteration with the inital guess given by forward
% Euler.

% Issues:
% Should use more accurate method for first step? - Trapezodial method or
% bdf1

% Open file to save data to (w=in write mode, +=wipe existing, t=text)
fileID = fopen('~/programming/dsmatlabtests/bdf2diagnostics','w+t');

% Create function to solve for BDF2 timestep
% bad??
bdf2_stepfun = @(t,yn,ynm1,ynp1) ( (2*h/3)*feval(odefun, t, ynp1) + (4/3)*yn - (1/3)*ynm1 );
bdf2_residual = @(t,yn,ynm1,ynp1) ( ynp1 - bdf2_stepfun(t,yn,ynm1,ynp1));

% Initialise
T_out = tspan(1):h:tspan(2);
N_eqns = length(y0);
Y_out = zeros(length(T_out),N_eqns);
Y_out(1,:) = y0;

% Calculate second value of y using forward Euler (bdf2 requires known
% values at 2 previous timesteps).
Y_out(2,:) = Y_out(1,:) + h*feval(odefun, T_out(2), Y_out(1,:));

for i= 3:length(T_out)
    % Solve the implicit timesteping equation using a rootfinder
    %Y_out(i,:) = simplerootfinder(bdf2_stepfun,h,T_out(i),Y_out(i-1,:),Y_out(i-2,:));
    Y_out(i,:) = newtonrootfinder(bdf2_residual,jacobian,h,T_out(i),Y_out(i-1,:),Y_out(i-2,:));
end


    function Y = simplerootfinder(stepfun,h,t,yn,ynm1)
        implicit_y_accuracy = 0.005;
        Y = yn + h*feval(odefun, t, yn); % use forward euler to get initial guess
        
        % Main while loop to solve for ynp1
        diff = 10*implicit_y_accuracy;   % ensure while loop is never accurate enough on first try
        j = 0;  % avoid infinite loops by counting number of iterations
        while (diff > implicit_y_accuracy)&&(j<10)
            Y_prev = Y;
            Y = feval(stepfun,t,yn,ynm1,Y);   % Find the next iteration
            diff = max( abs(Y(:) - Y_prev(:)) );    % Measure the change between steps.
            j = j+1;
        end
        
        % Give a warning if Y did not converge well enough
        if j>=10; warning('backward_euler:NCY',['The value of Y for t = ', ...
                num2str(t), ' did not converge to the required accuracy in 10 or less steps']); end
        
    end

    function Y = newtonrootfinder(residual,jacobian,h,t,yn,ynm1)
        % Solve the implicit function Y = residual(t,Y,yn,ynm1) for Y using
        % Newton's method.
        
        implicit_y_accuracy = 1e-5;
        Y = yn; % use previous step as initial guess
        
        % Main while loop to solve for new value of Y
        deltaY = 10*implicit_y_accuracy;
        k = 0;
        while (max(abs(deltaY)) > implicit_y_accuracy)&&(k<10)
            r = feval(residual,t,yn,ynm1,Y)';   %Evaluate residual and pass out H_eff
            J = feval(jacobian,t,Y);      % Evaluate jacobian with H_eff
            
            % Output data to file of fileID
            fprintf(fileID,'%d,%0.2e,%d,%0.2e\n',i,h,k,norm(r,2));
            
            deltaY = J\(-1*r);
            Y = Y + deltaY';
            k = k + 1;
        end
        
        if k>=10; warning('backward_euler:NCY',['The value of Y for t = ', ...
                num2str(t), ' did not converge to the required accuracy in 10 or less steps']); end
    end

end