function [i,j,k,f] = find_modes(room,f,c0)
    % This function help finding the closest frequency index 
    % for a certain mode considering the defined frequency range
    % it could be used as a tool with the green function for a room
    %
    % Parameters:
    % nx - specific modes in the x dimension
    % ny - specific modes in the y dimension
    % nz - specific modes in the z dimension
    % Optional parameters:
    %   room - dimension of the room [Lx,Ly,Lz] (default : [3.14, 4.38, 3.27] DTU room)
    %   f - frequency array, if used in addition to green function, 
    %       it should be the same for both, (default: 75 : 6.25e-2: (75+3200*6.25e-2) (same as meas))
    %     X-Axis first value:	 7.5000000000e+001	 
    %       X-Axis delta:	 6.2500000000e-002	 
    %
    % Returns:
    % idx - closest index to the specific mode accross the frequency axis
    % fm - frequency of the specific mode
    %
    % Exemple of usage of the function :
    %   room =[3.14, 4.38, 3.27];
    %   rS = [0,0,0];
    %   x = linspace(0,room(1),10);
    %   y = linspace(0,room(2),10);
    %   z = ones(size(x)) * room(3);
    %   rM = [x',y',z'];
    %   [G,f] = green_func_room(rM,rS,room, 'absorption', true);
    %   idx = find_f_modes(1,1,1);
    %   contourf(x,y,abs(G(:,:,idx)))
    arguments
        room = [3.28, 4.40, 3.28]
        f double = 700
        c0 double = 343
    end
    
    th = 0.1;
    nx_f = floor(2 *f * room(1) ./ (c0));ny_f = floor(2 *f * room(2) ./ (c0));nz_f = floor(2 *f * room(3) ./ (c0));
    nx = (nx_f-10):(nx_f+5);
    ny = (ny_f-10):(ny_f+5);
    nz = (nz_f-10):(nz_f+5);
    % Create meshgrid for all combinations
    [NX, NY, NZ] = meshgrid(nx, ny, nz);
    
    % Calculate frequencies for all combinations
    fm = c0 / 2 * sqrt((NX / room(1)).^2 + (NY / room(2)).^2 + (NZ / room(3)).^2);
    
    % Find indices within the frequency range
    valid_indices = find(fm >= f-50 & fm <= f+50);
    
    
    mode_idx = [NX(valid_indices), NY(valid_indices), NZ(valid_indices)];
    i = mode_idx(:,1);j = mode_idx(:,2);k = mode_idx(:,3);
    f = fm(valid_indices);

    % if fm > max(f)
    %     error("the frequency of the desired modes is greater than the frequency range defined, f="+round(fm,3)+" Hz")
    % end

    % disp("The frequency of the specific modes is : f="+round(fm,3)+" Hz")
    
    % [err, idx] = min(abs(f-fm));

    % if err > th
    %     disp("The frequency bins is "+round(err,4)+ "Hz away from the theorical frequency of the modes, consider defining a more precise frequency range")
    % end
end
