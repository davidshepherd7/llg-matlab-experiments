% Output script for the various LLG_solvers

% plot line showing path of M (parameterised by t)
figure
plot3(M_out(:,1),M_out(:,2),M_out(:,3),0,0,0,'o',M_out(1,1),M_out(1,2),M_out(1,3),'x')
legend('Path of magnetisation with time', 'Origin','Initial value','Location','NorthEast' )
xlabel('M_x'); ylabel('M_y'); zlabel('M_z');
movegui('northwest')

% also plot each component of M against t
figure
plot(T_out,M_out)
legend('Mx','My','Mz')
movegui('northeast')

% also plot each component of H against t
figure
hold on
plot(1*h:h:T_out(end),H_out)
legend('Hx','Hy','Hz')
movegui('southwest')