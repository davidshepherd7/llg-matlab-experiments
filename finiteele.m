function [y_approx, Ne] = finiteele(mesh_points, input_function, y0, y1)
%FINITEELE Calculate a finite element approximation for -y'' = F, y(0) =
%y0, y(1) = y1 over [0,1]. Uses Ne elements of the same size. F should be a
%function handle for a function with one variable.

% calculate Ne and h's from mesh points vector
Ne = length(mesh_points) - 1;
h = mesh_points(2:end) - mesh_points(1:end-1);


% initialising:
A_elements = cell(Ne,1);
f_elements = cell(Ne,1);
A_global = zeros(Ne+1,Ne+1);
f_global = zeros(Ne+1,1);


for e = 1:Ne    % for each element
    A_elements{e} = zeros(2);   %initialising
    f_elements{e} = zeros(2,1);
    
    % calculate matrix and f vector
    for i = 1:2
        % create function to be integrated in calculation of f_elements (phi_i * input_function)
        
        % error is probably coming from this part:
        func_to_integrate = @(x) (  ( x - mesh_points(e + (2-i) ) ) .* input_function(x)    );
        f_elements{e}(i) = (1/(h(e))) .* ((-1)^i) .* quad(func_to_integrate,mesh_points(e),mesh_points(e+1));
        
    end
    
    A_elements{e} = [1/h(e), -1/h(e); -1/h(e), 1/h(e)];
    
    
    % map local element matrix to global matrix
    A_global(e:e+1,e:e+1) = A_global(e:e+1,e:e+1) + A_elements{e};
    f_global(e:e+1) = f_global(e:e+1) + f_elements{e};
end


% apply boundary conditions
f_global(1) = y0;
A_global(1,:) = [1,zeros(1,Ne)];

f_global(end) = y1;
A_global(end,:) = [zeros(1,Ne),1];

% use sparse matrix since A_global is mostly zeros
A_global = sparse(A_global);

y_approx = A_global\f_global;  % approximation for y

f_finiteele = 1/h(1) * f_global

