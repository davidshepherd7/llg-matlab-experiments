function [dMdt] = LLG(M,H,Ms,gamma,alpha)
% The Landau-Lifshitz-Gilbert equation

% Set default values for constants:
if nargin < 3; Ms = 1; end
if nargin < 4; gamma = 0.2; end
if nargin < 5; alpha = 0.7; end

% Calculate dMdt
dMdt = (gamma/(1+alpha^2))*(cross(M,H) - (alpha/Ms)*cross(M,cross(M,H)));