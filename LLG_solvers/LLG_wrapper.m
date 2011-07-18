% Wrapper for LLG solvers to convert SI to arbitary units

% Saturation magnetisation:
Ms_SI = Ms;
Ms = 1;


% Applied field:
H_applied = H_applied/Ms_SI;
if length(H_applied) <= 3
    H_applied = repmat(H_applied,N,1);   % If constant convert to constant vector
end


% K1 conversion s.t. K1 = H_k ~ 1:
K1_SI = K1;
mu0 = 4*pi*1e-7;
K1 = K1*mu0;    % Re-normalise to give K1 ~ 1
K1 = 2*K1;  % Now H_k = K1


% Run function
[T_out,M_out,H_out] =  LLG_solver_simple(N,h,M0,H_applied,Ms,K1,easyaxis_direction);



% Convert back to original unit system
% To have gamma of order 0.1 we use time units of microseconds
T_out = T_out*1e-6;
h = h*1e-6;

% Just mulitply by Ms to convert back
M_out = M_out*Ms_SI;
H_out = H_out*Ms_SI;