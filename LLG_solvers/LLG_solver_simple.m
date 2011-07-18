function [T_out,M_out,H_out] =  LLG_solver_simple(N,h,M0,H_applied,Ms,K1,easyaxis_direction)
% Solve the LLG function for the simple case of a single, uniformly
% magnetised ellipsoid of rotation with uniaxial crystalline anisotropy.

% N is the number of times that the demagnetising field is recalculated. h
% is the step size between field recalculations. Note that this is not the
% same as the step size used for solution of the ODE which is decided by
% the built in matlab function ode45.

% H_applied should be a (3xN) matrix containing a list of the values of the
% applied field at each timestep.

% The magnetisation of the ellipsoid is assumed to always be saturated at
% Ms (to help with this M0 is normalised to length Ms before starting).

% A tiny 'kick' field is added to ensure the magnetisation cannot sit at 
% an unstable fixed point.

%==========================================================================
% Inputs
%==========================================================================

% Constants:
% γ = 2.21e5 m/As from: Application of the stereographic projection to studies of magnetization dynamics described by the Landau–Lifshitz–Gilbert equation, Journal of Physics A: Mathematical and Theoretical
gamma = 0.221; % in (m/A)*(1/micro sec)
alpha = 0.7;    % The relative strength of damping is propotional to alpha.

% Geometrical properties:
ellipsoid_axis_a = 2;   % The lengths of the axes of the ellipsoid.
ellipsoid_axis_b = 1;

% Kick field:
% Give a tiny 'kick' in all directions to prevent the magnetisation from
% sitting at a maximum energy. [This could be made more analogous to thermal
% effects by randomising the direction at each field recalculation step.]
H_kick = [1 1 1].*1e-6;

% Ensure unit vectors where applicable:
M0 = unit_vec(M0)*Ms;
easyaxis_direction = unit_vec(easyaxis_direction);


%==========================================================================
% Calculations
%==========================================================================

% Initialise output vectors:
H_out = zeros(N,3); 
M_out = M0;
T_out = 0;

for i = 1:N
    % Field (re)calculations:
    % Demag field:
    H_demag = ellipsoid_demag(ellipsoid_axis_a, ellipsoid_axis_b, M_out(end,:));
    % Find correct direction for crystalline anisotropy field using sign of dot
    % product then multiply by strength of field. [because of choice of units K1 = |H_k|]
    H_cryst = K1 * (sign(dot(M_out(end,:),easyaxis_direction)) * easyaxis_direction);
    % Add up total field:
    H_out(i,:) = H_demag + H_cryst + H_applied(i,:) + H_kick;
    H = H_out(i,:);
    
    % (Re)define LLG function with the new value of H:
    % [transpose taken because ode45 requires a column vector output]
    %dMdt = @(t,M) (gamma/(1+alpha^2))*(cross(M,H(i,:))...
    %   - (alpha/Ms)*cross(M,cross(M,H(i,:))))';
    % [cross product is main bottleneck of program so done explicitly for increased speed]
    dMdt = @(t,M) (gamma/(1+alpha^2))*[ ...
          M(2)*H(3) - M(3)*H(2) - (alpha/Ms)*( M(2)*(M(1)*H(2) - M(2)*H(1)) - M(3)*(M(3)*H(1) - M(1)*H(3))); ...  % dM_x/dt
          M(3)*H(1) - M(1)*H(3) - (alpha/Ms)*( M(3)*(M(2)*H(3) - M(3)*H(2)) - M(1)*(M(1)*H(2) - M(2)*H(1))); ...  % dM_y/dt
          M(1)*H(2) - M(2)*H(1) - (alpha/Ms)*( M(1)*(M(3)*H(1) - M(1)*H(3)) - M(2)*(M(2)*H(3) - M(3)*H(2)))];     % dM_z/dt

    % Run the ode solver for timestep from h*(i-1) to h*i:
    [T_solver,M_solver] = ode45(dMdt,[h*(i-1) h*i], M_out(end,:));
    
    % Normalise M vectors to have length Ms
    % [Probably can be commented out for smallish number of iterations]
    M_solver = Ms.*unit_vec(M_solver);
    
    % Store all T,M values from ode45:
    T_out = [T_out;T_solver];
    M_out = [M_out; M_solver];  
end