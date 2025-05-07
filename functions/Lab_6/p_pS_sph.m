function p = p_pS_sph(m, r, the, k, a, varargin)
    %   Calculates pressure field from a point source placed on a rigid sphere
    %    of radius a with modal order m at positions given 
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
    %       'norm_ax' - compute the normalized pressure relatvie to the  on-axis pressure (default: false)
    %       'amp'     - calculate actual amplitude (default: true)
    %
    %   Output:
    %       p       - complex pressure field [Pa]
    p = inputParser;
    addRequired(p, 'm');
    addRequired(p, 'r');
    addRequired(p, 'the');
    addRequired(p, 'k');
    addRequired(p, 'a');
    addParameter(p, 'Q', 1e-3);
    addParameter(p, 'norm_ax', false);
    addParameter(p, 'amp', true);
    
    parse(p, m, r, the, k, a, varargin{:});
    Q = p.Results.Q;
    norm_ax = p.Results.norm_ax;
    amp = p.Results.amp;

    c = 343;          
    rho = 1.21; 
    
    cth = cos(the);
    kr = k*r.';
    omega = k * c;
    [~,idxthe] = min(abs(the));
    Uthe = zeros(1,length(the));
    Uthe(idxthe) =1;
    p = zeros(length(r),length(the));
    for ii = 1:length(m)
        disp(m(ii))
        P = legendre(m(ii),cth);
        Pm = P(1,:);

        if amp == true
            Um = (m(ii)+1/2) * Q/ (2*pi*a^2) * Uthe; % not totally sure about this one
            Am = -1j * omega * rho *  Um / dr_sphankel2(m(ii),k*a);
        else
            Am = 1;
        end


        p = p +  Am .* Pm .*sphankel2(m(ii),kr); % time dependencies is not implemented on purpose
        
    end

end