function G = Gf_duct(duct_coord,rS,zM,duct, f,N_modes,mode_list, method)
    % Compute the Green's function for a rectangular duct.
    % Uses modal decomposition of a 3D plane wave.

    % Parameters:
    % duct_coord : [x, y] Coordinates of the observation point inside the duct.
    % rS : [x, y, z] Coordinates of the source point.
    % zM : array of z-coordinates where the Green's function is computed.
    % duct : [lx, ly] Dimensions of the rectangular duct.
    % f : frequency or array of frequencies of interest.
    % N_modes : [N, M] Maximum number of modes in x and y directions.
    % mode_list : [n, m] Specific modes to compute (optional).
    % method : Integer (1 or 2), determines the method of computation.
    
    % Returns:
    % G : Green's function values for the given frequencies and positions.
    
    arguments
        duct_coord 
        rS
        zM
        duct
        f
        N_modes = [20,20]; % Default mode numbers if not specified.
        mode_list = []; % Default: compute all modes.
        method = 1; % Default method.
    end

    % Constants
    c0 = 343;  % Speed of sound in air (m/s)
    omega = 2*pi *f;  % Angular frequency
    k = omega / c0;  % Wavenumber

    % Constrain the number of modes based on frequency
    N_modes = constrain_modes(N_modes,duct,f);
    % Compute the duct cross-sectional area
    S = prod(duct);
    lx = duct(1); 
    ly = duct(2);

    % Check if duct_coord contains 2D coordinates (x, y)
    if length(duct_coord) == 2
        if method == 1
            % Method 1: Compute Green's function using modal decomposition
            
            % Generate mode indices (nx, ny) based on N_modes
            if isempty(mode_list)
                [nx, ny] = ndgrid(0:N_modes(1), 0:N_modes(2));
                nx = nx(:);
                ny = ny(:);
                
                % Mode coefficients for normalization
                en = ones(size(nx));
                en(2:end) = en(2:end) .* 2;
                em = ones(size(ny));
                em(2:end) = em(2:end) .* 2;
                
                % Compute wave numbers in x and y directions
                coefx = nx*pi/lx; 
                coefy = ny*pi/ly;
            
                % Compute mode shapes at source and receiver positions
                phim_rS = sqrt(em.*en) .* cos(coefx .* rS(1)) .* cos(coefy.* rS(2));
                phim_rM = sqrt(em.*en) .* cos(coefx.* duct_coord(1)) .* cos(coefy.* duct_coord(2));
            else
                % Use specific mode list if provided
                mode_list = mode_list';
                if size(mode_list,2) == 2
                    [nx, ny] = ndgrid(mode_list(:,1),mode_list(:,2));
                else
                    mode_list = cat(2, mode_list, zeros(size(mode_list)));
                    [nx, ny] = ndgrid(mode_list(:,1), mode_list(:,2));
                end
                
                nx = nx(:);
                ny = ny(:);
                en = ones(size(nx)) * 2;
                em = ones(size(ny)) * 2;

                % Correct normalization for fundamental modes
                if mode_list(1,1) == 0
                    en(1) = 1;
                end
                if mode_list(1,2) == 0
                    em(1) = 1;
                end
                
                coefx = nx*pi/lx; 
                coefy = ny*pi/ly;
                    
                % Compute mode shapes at source and receiver positions
                phim_rS = sqrt(em.*en) .* cos(coefx .* rS(1)) .* cos(coefy.* rS(2));
                phim_rM = sqrt(em.*en) .* cos(coefx.* duct_coord(1)) .* cos(coefy.* duct_coord(2));
            end

            % Initialize Green's function array
            G = zeros(length(zM), length(f));
            
            % Compute Green's function for each frequency
            for ii = 1:length(f)
                kzmn = -sqrt(k(ii).^2 - coefx.^2 - coefy.^2); 
                G(:,ii) = sum((phim_rS .* phim_rM ./ kzmn) .* exp(-1j * kzmn .* (zM - rS(3))));
            end

        else
            % Method 2: Alternative approach to compute Green's function
            G = zeros(length(zM), length(f));
            for nn = 0:N_modes(1)
                for mm = 0:N_modes(2)
                    em = 2 - (mm == 0);
                    en = 2 - (nn == 0);
                    
                    coefx = nn * pi / lx; 
                    coefy = mm * pi / ly;
                    
                    phim_rS = sqrt(em .* en) .* cos(coefx .* rS(1)) .* cos(coefy .* rS(2));
                    phim_rM = sqrt(em .* en) .* cos(coefx .* duct_coord(1)) .* cos(coefy .* duct_coord(2));

                    kzmn = -sqrt(k.^2 - coefx.^2 - coefy.^2); 

                    for ii = 1:length(f)
                        G(:,ii) = G(:,ii) + ((phim_rS .* phim_rM ./ kzmn(ii)) .* exp(-1j * kzmn(ii) .* (zM' - rS(3))));
                    end
                end
            end
        end

    else
        % If a full cross-section is given, enforce scalar zM
        assert(isscalar(zM), "You need to fix a point in the tube to observe pressure across cross section")

        % Compute mode indices if mode_list is not provided
        if isempty(mode_list)
            [nx, ny] = ndgrid(0:N_modes(1), 0:N_modes(2));
            nx = nx(:);
            ny = ny(:);
            
            en = ones(size(nx));
            en(2:end) = en(2:end) .* 2;
            em = ones(size(ny));
            em(2:end) = em(2:end) .* 2;
            
            coefx = nx * pi / lx; 
            coefy = ny * pi / ly;
        end

        % Compute Green's function for all positions in the cross-section
        phim_rS = sqrt(em .* en) .* cos(coefx .* rS(1)) .* cos(coefy .* rS(2));
        G = zeros(length(duct_coord(:,1)), length(duct_coord(:,2)), length(f));
        
        for ll = 1:length(duct_coord(:,1))
            for jj = 1:length(duct_coord(:,2))
                phim_rM = sqrt(em .* en) .* cos(coefx .* duct_coord(ll,1)) .* cos(coefy .* duct_coord(jj,2));

                for ii = 1:length(f)
                    kzmn = -sqrt(k(ii).^2 - coefx.^2 - coefy.^2);
                    G(ll,jj,ii) = sum((phim_rS .* phim_rM ./ (kzmn+1e-9)) .* exp(-1j * kzmn .* (zM - rS(3))));
                end
            end
            disp("Processing Green's function for cross-section: " + ll + "/" + length(duct_coord(:,1)))
        end
    end

    % Normalize Green's function
    G = -1j / S * G;             
end

% Function to limit the number of modes based on frequency
function N_modes = constrain_modes(N_modes, duct, f)
    if isscalar(f)
        fmax = f;
    else
        fmax = f(end);
    end
       
    c0 = 343;
    nx_max = fmax/c0 * 2*pi * duct(1);
    ny_max = fmax/c0 * 2*pi * duct(2);

    f_mode_max = c0/2 * sqrt((nx_max/duct(1))^2 + (nx_max/duct(2))^2);
    
    while f_mode_max > fmax
        nx_max = nx_max - 1;
        ny_max = ny_max - 1;
        f_mode_max = c0/2 * sqrt((nx_max/duct(1))^2 + (nx_max/duct(2))^2);
    end

    N_modes = min(N_modes, [nx_max, ny_max]);
end
