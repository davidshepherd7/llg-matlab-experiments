function [ H_cryst_ineasydirection] = crystalline_anisotropy(K1, Ms)
% Calculate the effective field for a uniaxial crystalline anisotropy with
% K1, M, Ms as provided.

% using arbitary units: K1 = K1_true*mu0 so we set mu0 to 1
mu0 = 1;

H_cryst_ineasydirection = (2*K1)/(mu0*Ms);     % From Craik, pg 91 (units: A/m)
end

