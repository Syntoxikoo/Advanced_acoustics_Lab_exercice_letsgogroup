function Ynm = sph_harmonic(l, m, theta, phi)
    %   Calculates spherical harmonic Y_l^m for degree l and order m
    %   at spherical coordinates (theta, phi).
    %
    %   Inputs:
    %       l       - degree (non-negative integer)
    %       m       - order (integer, -l <= m <= l)
    %       theta   - polar (inclination) angle [rad] (array)
    %       phi     - azimuthal angle [rad] (array)
    %
    %   Output:
    %       Ynm     - spherical harmonic values
    
    
    if l < m 
        error("Invalid parameters: l must be ≥ m");
    end
    
    
    fact = factorial(l - m) / factorial(l + m);
    C = sqrt((2*l + 1)/(4*pi) * fact);
    
    Plm = legendre(l, cos(phi));
    if l ~= 0
        Plm = squeeze(Plm(m+1, :, :));   
    end
    
    
    Ynm = C * Plm .* exp(1j * m * theta);
end
