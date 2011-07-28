function magavi(T_out,M_out,start,finish,filename)
% Create video of the magnetisation changes (from M_out) and save as an avi
% file. Entries of M_out from start to finish are used in the video. An avi
% file with name filename_compressed is saved in the current directory
% (does not need to include .avi extension).

% Problem: video speed depends on the time steps used by the solver.
% ?? Call LLG_solver... from here and specify the times when we want M_out
% to the ode solver then find those times in the matrix and use them.

% Problem: have to manually find the "interesting" part.

clear video_frames
figure

% For each frame used:
for i = start : finish 
    % plot this frame
    clf
    plot3(M_out(:,1),M_out(:,2),M_out(:,3),0,0,0,'o',M_out(1,1),M_out(1,2),M_out(1,3),'x')
    arrow([0 0 0],[M_out(i,1),M_out(i,2),M_out(i,3)])
    %axis([-1 1 -1 1 -1 1]);
    xlabel('M_x'); ylabel('M_y'); zlabel('M_z');
    
    title(['t = ', num2str(T_out(i),3)]); % put the time in title

    % and add to list of frames.
    video_frames(i - (start -1)) = getframe(gcf);
end

% Convert to avi and save.
movie2avi(video_frames,filename);

% Compress avi file using ffmpeg (system command)
compress_command = ['!ffmpeg -i "', filename, '.avi" "', filename,'_compressed.avi"'];
eval(compress_command);

% Remove original (very large) file (system command)
rm_command = ['!rm "', filename, '.avi"'];
eval(rm_command);
end