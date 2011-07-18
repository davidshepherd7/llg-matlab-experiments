function [T_out,M_out,H_cart] =  LLG_solver_sph(N,h,M0,H_applied,Ms,K1,easyaxis_direction)
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

% A tiny random direction 'kick' field is added to ensure the magnetisation
% cannot sit at an unstable fixed point. To make more realistic should
% check temperature dependence.

%==========================================================================
% Inputs
%==========================================================================

% Choice of spherical polar system (for conversions)
coordsystem = 'elevation';

% Constants:
% gamma is the electron gyromagnetic ratio 
% = 2.21e5 m/As Application of the stereographic projection.., JoP A
gamma = 0.221; % in (m/A)*(1/micro sec)
alpha = 0.7;    % The relative strength of damping is propotional to alpha.

% Geometrical properties:
ellipsoid_axis_a = 2;   % The lengths of the axes of the ellipsoid.
ellipsoid_axis_b = 1;

% Convert M0 to spherical polars and ensure length = Ms:
M0 = carttosph(M0,coordsystem); M0(1) = Ms;

%==========================================================================
% Calculations
%==========================================================================

% Initialise output vectors:
% M_sph and T_out store data for all times from ode45.
% M_cart, H_cart, H_sph only store data from between runs of ode45.
H_cart = zeros(N,3); H_sph = H_cart; M_cart = H_cart;
M_cart(1,:) = sphtocart(M0,coordsystem);
M_sph = M0(2:3);
T_out = 0;

for i = 1:N
    % Field (re)calculations:
    % Demag field:
    H_demag = ellipsoid_demag(ellipsoid_axis_a, ellipsoid_axis_b, M_cart(end,:));
    % Find correct direction for crystalline anisotropy field using sign of dot
    % product then multiply by strength of field. [because of choice of units K1 = |H_k|]
    H_cryst = K1 *(sign(dot(M_cart(end,:),easyaxis_direction)) * easyaxis_direction);
    H_kick = (1e-3)*rand(1,3);
    % Add up total field:
    H_cart(i,:) = H_demag + H_cryst + H_applied + H_kick;
    
    % Convert field to spherical polars ready for time stepping
    H_sph(i,:) = carttosph(H_cart(i,:),coordsystem);
    
    % (Re)define LLG function with the new value of H:
    % For this section (only) M(1) = theta, M(2) = phi since no
    % need for r component.
    dMdt = @(t,M) (-gamma/(1+alpha^2)) *[ ...   %pre-factor
        ... %d(theta)/dt component:
        alpha*( H_cart(i,1)*cos(M(1))*cos(M(2)) + H_cart(i,2)*cos(M(1))*sin(M(2)) - H_cart(i,3)*sin(M(1))) ... % alpha* (H in theta direction)
        + (cos(M(2)*H_cart(i,2)) -sin(M(2))*H_cart(i,1))  ... % H in phi direction
        ; ...
        ... % d(phi)/dt component:
        (1/limited_sin(M(1))) * ( ...         % pre-factor for this component only
        - ( H_cart(i,1)*cos(M(1))*cos(M(2)) + H_cart(i,2)*cos(M(1))*sin(M(2)) - H_cart(i,3)*sin(M(1))) ...  % minus (H in theta direction)
        + alpha*(cos(M(2)*H_cart(i,2)) -sin(M(2))*H_cart(i,1)) ... % + alpha* H in phi direction
        )];

    % Run the ode solver for timestep from h*(i-1) to h*i:
    [T_solver,M_solver] = ode45(dMdt,[h*(i-1) h*i], M_sph(i,:));

    % Store results from ode solver for output and next loop:
    M_sph = [M_sph; M_solver];
    T_out = [T_out; T_solver];
    
    % Convert most recent M to cartesian coordinates ready for field calculations
    M_cart(i+1,:) = sphtocart([Ms, M_sph(end,:)],coordsystem);
    
end

% Convert M_sph to cartesian co ordinates for consistency with other solver
M_out = zeros(length(M_sph),3);
for i=1:length(M_sph)
    M_out(i,:) = sphtocart([Ms,M_sph(i,:)],coordsystem);
end