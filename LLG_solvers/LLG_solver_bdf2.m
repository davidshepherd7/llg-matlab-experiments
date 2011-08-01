function [T_out,M_out,H_out] =  LLG_solver_bdf2(T,M0,H_applied,Ms,K1,easyaxis_direction,alpha,h)
% Solve the LLG function for the simple case of a single, uniformly
% magnetised ellipsoid of rotation with uniaxial crystalline anisotropy.

% The magnetisation of the ellipsoid is assumed to always be saturated at
% Ms (to help with this M0 is normalised to length Ms before starting).

% A tiny random direction 'thermal' field is added to ensure the magnetisation
% cannot sit at an unstable fixed point. To make more realistic should
% check temperature dependence.

% List of differences from LLG_solver_simple:
% - Have to specify h (for now) since my solver is not adaptive.
% - Can work entirely in row vectors; ode45 only accepts column vectors but
% outputs rows!

%==========================================================================
% Inputs
%==========================================================================

% Need effective field as global since needed all over the place (still not
% happy with this)
global H;

% Constants:
% gamma is the electron gyromagnetic ratio 
% = 2.21e5 m/As - Application of the stereographic projection.., JoP A
gamma = 0.221; % in (m/A)*(1/micro sec)

% Geometrical properties:
ellipsoid_axis_a = 2;   % The lengths of the axes of the ellipsoid.
ellipsoid_axis_b = 1;

% Ensure vectors are unit vectors where applicable:
M0 = unit_vec(M0)*Ms;
easyaxis_direction = unit_vec(easyaxis_direction);

% Load jacobian function from file and subs in values of constants

jacobian = @(t,M) jacobian_bdf2_LLG(H,Ms,M,alpha,gamma);

%==========================================================================
% Nested functions
%==========================================================================
% function to calculate effective field
    function H_eff = H_eff_calc(t,M)
        % Find correct direction for crystalline anisotropy field using sign of dot
        % product then multiply by strength of field. [because of choice of units K1 = |H_k|]
        H_cryst = K1 * (sign(dot(M,easyaxis_direction)) * easyaxis_direction);
        
        H_demag = ellipsoid_demag(ellipsoid_axis_a, ellipsoid_axis_b, M); % Demag field
        H_thermal = (1e-3)*rand(1,3);        % Thermal field
        H_eff = H_demag + H_cryst + H_applied(t) + H_thermal;   % Add up total field:
        
    end

% function to compute dMdt
    function [dMdt H] = LLG_eqn_nested(t,M)
        % First calculate H-fields:
        H = H_eff_calc(t,M');
        
        % Now calculate dMdt using LLG function with the value of H just calculated:
        % [cross product is main bottleneck of program so done explicitly
        % for increased speed]
        dMdt = (gamma/(Ms*(1+alpha^2)))*[ ...
            M(2)*H(3) - M(3)*H(2) - (alpha/Ms)*( M(2)*(M(1)*H(2) - M(2)*H(1)) - M(3)*(M(3)*H(1) - M(1)*H(3))), ...  % dM_x/dt
            M(3)*H(1) - M(1)*H(3) - (alpha/Ms)*( M(3)*(M(2)*H(3) - M(3)*H(2)) - M(1)*(M(1)*H(2) - M(2)*H(1))), ...  % dM_y/dt
            M(1)*H(2) - M(2)*H(1) - (alpha/Ms)*( M(1)*(M(3)*H(1) - M(1)*H(3)) - M(2)*(M(2)*H(3) - M(3)*H(2)))];     % dM_z/dt
    end

%==========================================================================
% Calculations
%==========================================================================

% Run the ode solver
[T_out,M_out] = bdf2(@LLG_eqn_nested,[0 T], M0, h, jacobian);

% ??: M vectors are NOT normalised, not sure how to implement this.

% If H_out is requested calculate H values for each (M,t) from ode solver.
% [?? inefficient but other ways to do this are dodgy or require modifying
% the ode solver, vectorising is very messy so just use loop]
if nargout >2
    H_out = zeros(length(T_out),3);
    for i=1:length(T_out)
        H_out(i,:) = H_eff_calc(T_out(i), M_out(i,:));
    end
end

end

