function [T_out,M_out,H_out] =  LLG_solver_simple(N,h,M0,H_applied,Ms,K1,easyaxis_direction)
% Solve the LLG function for the simple case of a single, uniformly
% magnetised ellipsoid of rotation with uniaxial crystalline anisotropy.

% N is the number of times that the demagnetising field is recalculated. h
% is the step size between field recalculations. Note that this is not the
% same as the step size used for solution of the ODE which is decided by
% the built in matlab function ode45.

% The magnetisation of the ellipsoid is assumed to always be saturated at
% Ms (to help with this M0 is automatically normalised to have length Ms).

% The program adds in a tiny 'kick' field to ensure the magnetisation
% cannot sit at an unstable fixed point.

% Using the SI unit system:
% Units of H_applied are A/m
% Units of a and b are irrelevant (so long as they are the same)
% Ms are in A/m
% K1 are in J/m^3
% time units (e.g. for h) are femtoseconds

%??????????????????????????????????????????????????????????????????????????
% Anisotropy needs to go both directions (both types)
% Finish checking units
% Ellipsoid demag calcualations wrong?
%??????????????????????????????????????????????????????????????????????????


%==========================================================================
% Inputs
%==========================================================================

% Constants:
gamma = 0.17e-4 ;% The gyromagnetic ratio - test value that works, equivalent to working in femtoseconds
%gamma = 1.760859e11; % real value of gamma for electrons, too large for
%convergence?
alpha = 0.7;    % The relative strength of damping is propotional to alpha.
mu0 = 4*pi*1e-7;    % Magnetic permeabillity of vacuum (units: Vs / A(m^2))

% Geometrical properties:
ellipsoid_axis_a = 2;   % The lengths of the axes of the ellipsoid.
ellipsoid_axis_b = 1;

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

% Calculate crystalline anisotropy effective field
H_cryst_easyaxis = (2*K1)/(mu0*Ms);      % From Craik, pg 91 (units: A/m)
H_cryst = H_cryst_easyaxis*easyaxis_direction;

% Initial values for output vectors
% [Note that M_out and T_out are not of a pre-defined size because their
% size is dependent on the number of steps used by ode45, which is hard to
% predict.]
M_out = M0;
T_out = 0;
H_out = H_applied;

for i = 1:N
    
    % Demag and total field (re)calculation
    % Units for both: A/m
    H_demag = ellipsoid_demag(ellipsoid_axis_a, ellipsoid_axis_b, M_out(end,:));
    H = H_demag + H_cryst + H_applied + H_kick;
    
    H_out = [H_out; H];
    
    % (Re)define LLG function with the new value of H
    % [transpose taken because ode45 requires a column vector output]
    dMdt = @(t,M) (gamma/(1+alpha^2))*(cross(M,H)...
        - (alpha/Ms)*cross(M,cross(M,H)))';

    % Run the ode solver for timestep from h*(i-1) to h*i
    [T_solver,M_solver] = ode45(dMdt,[h*(i-1) h*i], M_out(end,:));

    % Store results from ode solver for output and next loop
    M_out = [M_out; M_solver];
    T_out = [T_out; T_solver];
end

H_out = H_out(2:end,:);