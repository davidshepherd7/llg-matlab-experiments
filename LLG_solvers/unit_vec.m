function x_unit = unit_vec(x)
% Normalise a list of vectors (in the form below) to unit vectors
% [ x1 x2 x3 ... ;
%   y1 y2 y3 ... ;
%   ........ ... ;
%   z1 z2 z3 ... ]

% Initialise output
x_unit = zeros(size(x));

% For each vector
for i=1:size(x(:,1))
    
    % Calculate normlising factor
    d = 1/norm(x(i,:),2);
    
    % Check if close to dividing by zero
    if 1/d < 1e-20;
        warning('unit_vec:small_norm', 'The norm of the vector is very small, this may be close to division by zero')
    end
    
    % Normalise to unit vector
    x_unit(i,:) = d*x(i,:);
end