function [H_loop_handle] = H_loop_gen(H_min,H_max,T,theta)
% Generate a loop of fields (in the z direction) from H_min to H_max and
% back again in time T. The field is at angle theta to the z
% axis.

    function H = H_loop(t)
        if t < T/2
            H_strength = H_min + t*(H_max - H_min)/(T/2);
        else
            H_strength = H_max - (t-T/2)*(H_max - H_min)/(T/2);
        end
        H = sphtocart([H_strength,theta,0]);
        
    end

H_loop_handle = @H_loop;
end