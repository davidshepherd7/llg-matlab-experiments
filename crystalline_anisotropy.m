function [ H_cryst] = crystalline_anisotropy(K1, easyaxis_direction, M, Ms, mu0 )
% Calculate the effective field for a uniaxial crystalline anisotropy with
% K1, M, Ms as provided.

H_cryst_easyaxis = (2*K1)/(mu0*Ms);      % From Craik, pg 91 (units: A/m)
H_cryst = H_cryst_easyaxis*easyaxis_direction;

% Ensure field points the right direction along the axis
if (M*easyaxis_direction' < 0)
    H_cryst = -1.*H_cryst;
end

end

