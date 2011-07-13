function [T_out,M_cart,H_cart] =  LLG_solver_sph(N,h,M0,H_applied,Ms,K1,easyaxis_direction)
% Solve the LLG function for the simple case of a single, uniformly
% magnetised ellipsoid of rotation with uniaxial crystalline anisotropy.

% N is the number of times that the demagnetising field is recalculated. h
% is the step size between field recalculations. Note that this is not the
% same as the step size used for solution of the ODE which is decided by
% the built in matlab function ode45.

% The magnetisation of the ellipsoid is assumed to always be saturated at
% Ms (to help with this M0 is normalised to length Ms before starting).

% A tiny 'kick' field is added to ensure the magnetisation cannot sit at 
% an unstable fixed point.

% Using the SI unit system:
% Units of H_applied are A/m
% Units of a and b are irrelevant (so long as they are the same)
% Ms is in A/m
% K1 is in J/m^3
% Units of time (e.g. for h, gamma) are nanoseconds


%==========================================================================
% Inputs
%==========================================================================

% Constants:
%gamma = 0.17e-4 ;% The gyromagnetic ratio - test value that works, equivalent to working in femtoseconds
%gamma = 1.760859e11 1/Ts; % real value of gamma for electrons, too large for convergence?
% γ = 2.21e5 m/As from: Application of the stereographic projection to studies of magnetization dynamics described by the Landau–Lifshitz–Gilbert equation, Journal of Physics A: Mathematical and Theoretical
gamma = 2.21e-4; % in (m/A)*(1/ns), probably more appropriate units than those using Tesla.

alpha = 0.7;    % The relative strength of damping is propotional to alpha.
mu0 = 4*pi*1e-7;    % Magnetic permeabillity of vacuum (units: Vs / A(m^2))

% Geometrical properties:
ellipsoid_axis_a = 2;   % The lengths of the axes of the ellipsoid.
ellipsoid_axis_b = 1;

% Kick field:
% Give a tiny 'kick' in all directions to prevent the magnetisation from
% sitting at a maximum energy. [This could be made more analogous to thermal
% effects by randomising the direction at each field recalculation step.]
H_kick = [1 1 1].*1e-6;

% Convert M0 to spherical polars and ensure length = Ms:
M0 = sphtocart(M0); M0(1) = Ms;


%==========================================================================
% Calculations
%==========================================================================

% Initialise output vectors:
H_cart = zeros(N,3); M_sph = H_cart; H_sph = H_cart; M_cart = H_cart;
T_out = zeros(N,1);
M_sph(1,:) = M0;
M_cart(1,:) = carttosph([Ms,M_sph(1,:)]);

for i = 1:N
    
    % Field (re)calculations:
    % Units for all: A/m
    % [Recalculation of cryst is needed to ensure effective field H_cryst is
    % along the correct axis at all times.]
    H_demag = ellipsoid_demag(ellipsoid_axis_a, ellipsoid_axis_b, M_cart(end,:));
    H_cryst = crystalline_anisotropy(K1, easyaxis_direction, M_cart(end,:), Ms, mu0 );
    H_cart(i,:) = H_demag + H_cryst + H_applied + H_kick;
    
    % Convert field to spherical polars ready for time stepping
    H_sph(i,:) = carttosph(H_cart(i,:));
    
    % (Re)define LLG function with the new value of H:
    %????????????????? wrong!
    dMdt = @(t,M) (gamma/(1+alpha^2)) * [alpha*H_sph(i,1) + H_sph(i,2); ...
           (1/limited_sin(M_sph(i,1))) * ( H_sph(i,1) - alpha*H_sph(i,2) )];

    % Run the ode solver for timestep from h*(i-1) to h*i:
    [T_solver,M_solver] = ode45(dMdt,[h*(i-1) h*i], M_sph(end,2:3));

    % Store results from ode solver for output and next loop:
    M_sph(i+1,:) = [Ms M_solver(end,:)];
    T_out(i+1,:) = T_solver(end,:);

    % Convert M to cartesian coordinates ready for field calculations
    M_cart(i+1,:) = carttosph([Ms,M_sph(i,:)]);
end