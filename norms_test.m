
% test size of changes in length of M over time

Ms = zeros(length(M_out),1);
for i = 1:length(M_out)
    Ms(i) = norm(M_out(i,:),2);
end

plot(Ms)

% result: change is of order 1e-3 for M of order 1e6, insignificant