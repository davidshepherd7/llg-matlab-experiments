x_0 = 0.5
x_1 = 0.6

f = @(x) (x-x_0)/(x_1-x_0);

% 2 node approximation
% so x = +- srt(3)/3, w=1

x_eval = (sqrt(3)/3).*[-1,1];

x_eval = x_eval*((x_1-x_0)/2) + ((x_1 + x_0)/2);

plot([x_0:0.01:x_1],f([x_0:0.01:x_1]))

weights_2 = [1,1];

quadest =((x_1 - x_0)/2)*dot(f(x_eval),weights_2)