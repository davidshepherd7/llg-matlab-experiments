function [errors, N] = approximation_tests(input_function, exact_function,y0,y1)
% Test function finitediff and finiteele for various values of N. Output is
% the errors (with finite difference in the first column, finite element in the second)
% and the values of N. input_function and exact_function should be function
% handles.

j_max = 6;

errors = zeros(j_max,2);
N = zeros(j_max,1);

for j = 1:j_max;
    
    N(j) = 2^j +2;
    
    h = 1/(N(j)-1);
    mesh_points = 0:h:1;
    
    %calculate function values
    y_exact = exact_function(mesh_points)';
    
    y_finitediff = finitediff(input_function,N(j),y0,y1);
    y_finiteele = finiteele(mesh_points, input_function, y0, y1);
    
    %calculate errors
    errors(j,1) = max(abs(y_exact - y_finitediff));
    errors(j,2) = max(abs(y_exact - y_finiteele));
end


