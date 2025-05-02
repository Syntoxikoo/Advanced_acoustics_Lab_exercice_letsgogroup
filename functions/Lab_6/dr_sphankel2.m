function dh = dr_sphankel2(m, r)
    % Compute derivative of spherical Hankel function of the second kind
    % with respect to r
    if m ~= 0
        h_m1 = m * sphankel2(m-1,r);
    else
        h_m1 = 0;
    end
        dh = 1/(2*m+1) * (h_m1 - (m+1) * sphankel2(m+1,r));
    end