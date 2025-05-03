function [ur, uthe] = urad_sph_src(m, r, the, k, a, varargin)
    %   Calculates particle velocity from a spherical source of radius a 
    %   with modal order m at positions given 
    %   by spherical coordinates (r, the). The frequency of the wave is given in term
    %   of the wave number.
    %
    %   Inputs:
    %       m       - modal order (scalar or array)
    %       r       - radial coordinates (array)
    %       the     - polar angles [rad] (array)
    %       k       - wavenumber 
    %       a       - source radius
    %
    %   Optional parameters:
    %       Q         - volume velocity (default: 1e-3)
    %       'amp'     - calculate actual amplitude (default: true)
    %       'dthe"    - bool, if true compute uthe
    %
    %   Output:
    %       ur       - particle velocity [r,the] in the r directiong
    %       uthe       - particle velocity [r,the]  in the theta direction
    
    p = inputParser;
    addRequired(p, 'm');
    addRequired(p, 'r');
    addRequired(p, 'the');
    addRequired(p, 'k');
    addRequired(p, 'a');
    addParameter(p, 'Q', 1e-3);
    addParameter(p, 'amp', true);
    addParameter(p, 'dthe',false);

    
    parse(p, m, r, the, k, a, varargin{:});
    Q = p.Results.Q;
    amp = p.Results.amp;
    dthe = p.Results.dthe;

    c = 343;          
    rho = 1.21; 
    
    cth = cos(the);
    kr = k*r.';
    omega = k * c;

    ur = zeros(length(r),length(the));
    uthe = zeros(length(r),length(the));
    for ii = 1:length(m)
        P = legendre(m(ii),cth);
        Pm = P(1,:);

        if amp == true
            Um = (m(ii)+1/2) * Q/ (2*pi*a^2); % not totally sure about this one
            Am = -1j * omega * rho *  Um / dr_sphankel2(m(ii),a);
        else
            Am = 1;
        end
          
        ur = ur + (-Am) / (1j * rho * c) * Pm .*dr_sphankel2(m(ii),kr);

        if dthe == true
            dP = legendre_derivative(m(ii),the,P);
            uthe = uthe + (-Am) / (1j * rho * omega) * dP(1,:) .* (sphankel2(m(ii),kr)./ r);
    
        end
    end

end