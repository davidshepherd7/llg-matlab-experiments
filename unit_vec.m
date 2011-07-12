function x_unit = unit_vec(x)

if norm(x,2) < 1e-20;
    warning('unit_vec:small_norm', 'The norm of the vector is very small, this may be close to division by zero')
end

x_unit = x/norm(x,2);