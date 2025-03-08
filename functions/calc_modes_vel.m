function [Cp,Cg,nx] = calc_modes_vel(duct,f,N_modes)
    % duct : [lx, ly] Dimensions of the rectangular duct.
    arguments
    duct = [0.7, 1];
    f = 0:1000;
    N_modes = [8,8];
    end
    c0 = 343;
    lx = duct(1); 
    ly = duct(2);

    [nx, ny] = ndgrid(0:N_modes(1), 0:N_modes(2));
    nx = nx(:);
    ny = ny(:);
    

    % Compute wave numbers in x and y directions
    coefx = nx*pi/lx; 
    coefy = ny*pi/ly;

    k = 2*pi*f/c0;
    
    Cp = zeros(length(f),length(nx));
    Cg = zeros(length(f),length(nx));
    for ii = 1:length(f)
        kzmn = -sqrt(k(ii).^2 - coefx.^2 - coefy.^2);
        Cp(ii,:) = 2*pi*f(ii)./(kzmn+1e-6);
        Cg(ii,:) = c0 * sqrt(1 - (coefx.^2 + coefy.^2)./k(ii).^2);
    end
end