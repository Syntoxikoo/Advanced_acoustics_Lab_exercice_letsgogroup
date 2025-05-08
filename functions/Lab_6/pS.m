function p = pS(r,k,Q,a)
    % simulate the pressure radiated by a point source
    c = 343;
    rho = 1.21;

    omega = k * c;

    p = 1j *  omega * rho * (a^2*Q)./ ( (r-a)) .* exp(1j* k*(r-a));
    
end