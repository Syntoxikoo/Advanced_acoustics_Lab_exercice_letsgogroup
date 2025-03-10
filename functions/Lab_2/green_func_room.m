function [G,f] = green_func_room(rM, rS, room, varargin)
    % Compute the Green's function for a rectangular room
    % Uses modal decomposition of a 3D plane wave
    %
    % Parameters:
    % rM - Coordinates of the observation point [x,y,z]
    % rS - Coordinates of the source point [x,y,z]
    % room - Dimensions of the rectangular room [lx,ly,lz]
    % Optional parameters:
    %   'f' - Frequency array (default: 75 : 6.25e-2: (75+3200*6.25e-2) (same as meas))
    %   'absorption' - Boolean flag for absorption (default: false)
    %   'T60' - Reverberation time, could be either contstan or differents for each modes
    %   'compute' - number of modes computed 'range' or 'all' (default: 'range')
    %   'fc' - cutoff frequency of the mode for faster compute, default : 300
    %
    % Returns:
    % G - Green's function for given frequencies and positions 
    %     (if rM is a point in the room : [G] = (1 x Nfreq)), 
    %     (if rM is an array (cross section) : [G] = (Ndim1 x Ndim2 x Nfreq))
    %
    % Exemple of usage of the function for full frequency response :
    %   room =[3.14, 4.38, 3.27];
    %   rS = [0,0,0];
    %   x = room(1); y = room(2); z = 0
    %   rM = [x,y,z];
    %   [G,f] = green_func_room(rM,rS,room, 'absorption', true); 
    %   idx = find_f_modes(1,1,1);
    %   contourf(x,y,abs(G(:,:,idx)))
    %
    % Exemple of use case for cross section view :
    %   room =[3.14, 4.38, 3.27];
    %   rS = [0,0,0];
    %   x = linspace(0,room(1),10);
    %   y = linspace(0,room(2),10);
    %   z = ones(size(x)) * room(3);
    %   rM = [x',y',z'];
    %   [G,f] = green_func_room(rM,rS,room, 'absorption', true);
    %   idx = find_f_modes(1,1,1);
    %   contourf(x,y,abs(G(:,:,idx)))
    
    % -------------------------------- init function -----------------------------
    p = inputParser;
    addRequired(p, 'rM');
    addRequired(p, 'rS');
    addOptional(p, 'room', [3.14, 4.38, 3.27]);
    addOptional(p, 'f', 75 : 6.25e-2: (75+3200*6.25e-2));
    addParameter(p, 'absorption', false, @islogical);
    addParameter(p, 'T60', [3.5], @isnumeric);
    addParameter(p, 'compute', 'range');
    addParameter(p, "fc", 300, @isnumeric)
    
    parse(p, rM, rS, room, varargin{:});
    
    f = p.Results.f;
    absorption = p.Results.absorption;
    T60 = p.Results.T60;
    compute = (p.Results.compute);
    fc = p.Results.fc;

    % Constants
    c0 = 343;
    rho0 = 1.21;
    omega = 2 * pi * f;
    k = omega / c0;
    eps = 1e-9;

    lx = room(1); ly = room(2); lz = room(3);
    V = prod(room);
    % ------------------------------------------------------------------------------

    % Determine maximum mode in function of the frequency 
    % range or the defined cutoff
    if strcmp(compute, 'range')
        max_mode = compute_max_mode(fc, c0, room);
    else
        max_mode = compute_max_mode(f, c0, room);
    end
    
    %----------------------- Compute mode frequency --------------------------------
    [nx_grid, ny_grid, nz_grid] = ndgrid(0:max_mode, 0:max_mode, 0:max_mode);
    nx = nx_grid(:);
    ny = ny_grid(:);
    nz = nz_grid(:);
    
    pi_over_lx = pi / lx;
    pi_over_ly = pi / ly;
    pi_over_lz = pi / lz;
    
    fm = c0 / 2 * sqrt((nx / lx).^2 + (ny / ly).^2 + (nz / lz).^2);
    km = 2 * pi * fm / c0;

    % ------------------------------------------------------------------------------
    
    % add absorption with time constant, from specific reverb time 
    if absorption
        tau_m = T60/(6*log(10));
    else
        tau_m = 0;
    end  

    % --------------- Compute Gf for a paire of source - receiver ------------------
    if length(rM) <=3
        % Mode shapes at observation and source points 
        phim_rM = cos(nx * pi_over_lx * rM(1)) .* ...
                cos(ny * pi_over_ly * rM(2)) .* ...
                cos(nz * pi_over_lz * rM(3));
                
        phim_rS = cos(nx * pi_over_lx * rS(1)) .* ...
                cos(ny * pi_over_ly * rS(2)) .* ...
                cos(nz * pi_over_lz * rS(3));
                
        G = zeros(size(f));
        
  
        % compute green function 
        for ii = 1:length(f)
            G(ii) = sum(phim_rM .* phim_rS ./ (k(ii)^2 - km.^2 + eps - 1j*k(ii)/(tau_m*c0)));
        end
    % ------------------------------------------------------------------------------

    % --------------- Compute Gf for a the entire cross section --------------------
    else
        % ------ finding which dimension is constant -----
        uniqX = unique(rM(:,1));
        uniqY = unique(rM(:,2));
        uniqZ = unique(rM(:,3));
        nUniqPerDim = [length(uniqX), length(uniqY), length(uniqZ)];
        
        [~, constDim] = min(nUniqPerDim);
        cst = unique(rM(:,constDim));
        
        varyingDims = setdiff([1,2,3], constDim);
        dim1 = varyingDims(1);
        dim2 = varyingDims(2);
        r1 = unique(rM(:,dim1));
        r2 = unique(rM(:,dim2));

        [r1_g, r2_g] = ndgrid(r1,r2);
        rconst = ones(size(r1_g)) * cst;
        
        [nRows, nCols] = size(r1_g);
        numPoints = nRows * nCols;
            
        coords = zeros(numPoints, 3);
            
        switch constDim
            case 1  % X is constant
                coords(:,1) = rconst(:);
                coords(:,2) = r1_g(:);
                coords(:,3) = r2_g(:);
            case 2  % Y is constant
                coords(:,1) = r1_g(:);
                coords(:,2) = rconst(:);
                coords(:,3) = r2_g(:);
            case 3  % Z is constant
                coords(:,1) = r1_g(:);
                coords(:,2) = r2_g(:);
                coords(:,3) = rconst(:);
        end
        % ---------------------------------------------

        cos_x = cos(nx * pi_over_lx * coords(:,1)');
        cos_y = cos(ny * pi_over_ly * coords(:,2)');
        cos_z = cos(nz * pi_over_lz * coords(:,3)');
        phim_rM = cos_x .* cos_y .* cos_z;
        phim_rM = phim_rM.';
            
        phim_rS = cos(nx * pi_over_lx * rS(1)) .* ...
                cos(ny * pi_over_ly * rS(2)) .* ...
                cos(nz * pi_over_lz * rS(3));
                
        G = zeros(nRows, nCols, length(f));
        
        for ii = 1:length(f)
            if absorption
                den = (k(ii)^2 - km.^2 + eps - 1j*k(ii)/(tau_m*c0));
            else
                den = (k(ii)^2 - km.^2 + eps );
            end
            G_tmp = sum(phim_rM.* (phim_rS ./ den).',2);

            G(:,:,ii) = reshape(G_tmp, [nRows,nCols]);
        end
    end
    % ------------------------------------------------------------------------------
    G = G * (-1 / V);
end


function max_mode = compute_max_mode(f, c0, room)    
    if numel(f) > 1
        f_max = max(f);
    else
        f_max = f;
    end
    max_mode = ceil(f_max * 2 * min(room) / c0);
end