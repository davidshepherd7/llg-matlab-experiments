function [T_out,M_out,H_out] =  LLG_wrapper(T,M0,H_applied,Ms,K1,easyaxis_direction)
% Wrapper for LLG solvers to convert SI to arbitary units

% Saturation magnetisation:
Ms_SI = Ms;
Ms = 1;


% Applied field:
if isnumeric(H_applied)
    H_applied = @(t) (H_applied);   % If constant convert to constant fn
end

H_applied_normalised = @(t) (H_applied(t)/Ms_SI);


% K1 conversion s.t. K1 = H_k ~ 1:
K1_SI = K1;
mu0 = 4*pi*1e-7;
K1 = K1*mu0;    % Re-normalise to give K1 ~ 1
K1 = 2*K1;  % Now H_k = K1

% Set value of alpha (done here because easier to change ??bad place)
alpha = 0.7;


% Run function
[T_out,M_out,H_out] =  LLG_solver_simple(T,M0,H_applied_normalised,Ms,K1,easyaxis_direction,alpha);



% Convert back to original unit system
% ?? not sure on time units yet
T_out = T_out;

% Just mulitply by Ms to convert back
M_out = M_out*Ms_SI;
H_out = H_out*Ms_SI;