% Calculate the rotation of a single domain, magnetically soft 
%(no crystalline anisotropy), ellipsoid of rotation in an applied magnetic field.

% _cart = in cartesian co-ordinates
% _sph = in spherical polar co-ordinates

clear;

%-------------------------------------------------------------------------
% Set up coefficients and initial values
%-------------------------------------------------------------------------

M_s = 1; % saturation magnetisation
H_app_cart = [0,0,1]; % applied field in cartesian co-ordinates
M_0_sph = [0,0]; % initial magnetisation direction in spherical polars ??placeholder

% dimensions of ellipsoid of rotation
a = 2; 
b = 1;

% LLG parameters - defined in function LLG_sph until find how to define
% here.
% alpha = 0.8; % ??placeholder
% gamma = 0.5; % ??placeholder
% d = gamma/(1 + alpha^2); % used often so pre-calculate

% time step set up
t_step = 0.1; % time step to use 
N_time_steps = 10; % number of time steps to take


%-------------------------------------------------------------------------
% Calculate and convert fields
%-------------------------------------------------------------------------

H_demag_cart = [0, 0, ellipsoid_demag(a, b, M_s)] % calculate H_demag

H_eff_sph = cart2sph(H_demag_cart + H_app_cart) %convert to polars


%-------------------------------------------------------------------------
% Integrate LLG over time:
%-------------------------------------------------------------------------

% Needs a stiff solver (ode15s) when exhange coupling is included but for now
% non-stiff should be ok (ode45).

M_solution = ode45(@LLG_sph,[0, 10],M_0_sph);







