%Messing with magnetic energy equation

a = (-0.5:0.01:1.5);
theta = pi * a;

phi = 0;
K = 1;
V = 1;
H_k = 1;

H_eff = [0,0.5,1,1.2];
E = zeros(length(H_eff),length(a));

for i=1:length(H_eff)

E(i,:) = K*V*( (sin(theta)).^2  -  2* (H_eff(i)/H_k) .* cos(theta - phi) );

end

plot(a,E);
legend('H_{eff} = 0','0.5','1','1.2','Location','NorthWest')
xlabel('angle/ pi radians')
ylabel('Energy barrier/ some natural units')