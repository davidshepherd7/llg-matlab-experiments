% Output script for the various LLG_solvers

% plot line showing path of M (parameterised by t)
figure
plot3(M_out(:,1),M_out(:,2),M_out(:,3),0,0,0,'o',M_out(1,1),M_out(1,2),M_out(1,3),'x')
legend('Path of magnetisation with time', 'Origin','Initial value','Location','NorthEast' )
xlabel('M_x'); ylabel('M_y'); zlabel('M_z');

% also plot each component of M against t
figure
plot(T_out,M_out)
legend('Mx','My','Mz')
movegui('northeast')

% also plot each component of H against t
figure
plot(T_out,H_out)
legend('Hx','Hy','Hz')
movegui('southeast')

% plot a hystersis loop if H_applied is non-constant
if max(abs(H_applied(0) - H_applied(T/2 - 5))) > 0.0000001
    
    H_applied_values = zeros(length(T_out),3);  % calculate values of applied field
    for i=1:length(T_out)
        H_applied_values(i,:) = H_applied(T_out(i));
    end
    
    % plot z magnetisation against z applied field
    figure
    plot(H_applied_values(:,3),M_out(:,3))
    xlabel('z-applied field')
    ylabel('z-magnetisation')
    movegui('south')
    
end
    
