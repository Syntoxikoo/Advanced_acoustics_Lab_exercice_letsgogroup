function p = pS(r,k,Q)
    % simulate the pressure radiated by a point source
    c = 343;
    rho = 1.21;

    omega = k * c;

    p = 1j *  omega * rho * Q ./ (4 * pi * r) .* exp(1j* k*r);
    
end