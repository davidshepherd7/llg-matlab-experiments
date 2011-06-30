function [dM] = LLG_sph(M, H)
% The Landau-Lifschitz-Gilbert equation in spherical polar coordinates.
% M(1) = M_theta, M(2) = M_phi, similarly for H. M_r = |M_s| (fixed)

alpha = 0.8;
gamma = 0.5;

dM = (gamma/(1 + alpha^2)) * [ -alpha*H(1) - H(2); (1./sin(M(1))) .* (H(1) - alpha*H(2))];