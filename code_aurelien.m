% Parameters
c_0 = 343 ;
rho_0 = 1.2 ;

x_0 = 0 ;
y_0 = 0 ;
r_0 = [x_0, y_0, 0] ;
x = 0.5 ;
y = 0.5 ;
r = [x, y, 0] ;
a = 1 ; %height
b = 1 ; %width
S = a*b ; %cross-sectional area

n_max = 40 ;
m_max = 40 ;

%% fixed freq, over z

f = 10000 ;
omega = 2*pi*f ;
k = omega / c_0 ; %wavenumber 

z_max = 10;
z = linspace(0,z_max,1000) ;

G = zeros(1, length(z));

for iz = 1:length(z)
    sum = 0 ;
    z_val = z(iz);
    for m = 0:m_max
        for n = 0:n_max
            if m == 0
                eps_m = 1 ;
            else 
                eps_m = 2 ;
            end
            if n == 0
                eps_n = 1 ;
            else 
                eps_n = 2 ;
            end
            k_zmn = sqrt(k^2 - (m*pi/a)^2 - (n*pi/b)^2 ) ; 
            psi_mn_r = sqrt(eps_m*eps_n) * cos(m*pi*x/a) * cos(n*pi*y/a) ;
            psi_mn_r0 = sqrt(eps_m*eps_n) * cos(m*pi*x_0/a) * cos(n*pi*y_0/a) ;
            sum = sum + psi_mn_r.*psi_mn_r0./k_zmn .* exp(-1j.*k_zmn.*z_val) ;
        end
    end 
   G(iz) = -1j/S .* sum;
end

figure(1) ;
plot(z,20*log10((abs(G)/(2*10^-5)))) ;
grid on;
xlabel('Distance (m)') ;
ylabel('Pressure level (dB SPL)') ;
title('Pressure in function of distance') ;

%% fixed z, over freq.
f = linspace(1,1000,1000-1) ; %frequency (we take on point off 
% to not land on infinite points)

z = 0.1;

G = zeros(1, length(f));

for ifrq = 1:length(f)
    sum = 0 ;
    f_val = f(ifrq);
    omega = 2*pi*f_val ;
    k = omega / c_0 ;
    for m = 0:m_max
        for n = 0:n_max
            if m == 0
                eps_m = 1 ;
            else 
                eps_m = 2 ;
            end
            if n == 0
                eps_n = 1 ;
            else 
                eps_n = 2 ;
            end
            k_zmn = sqrt(k^2 - (m*pi/a)^2 - (n*pi/b)^2 ) ; 
            psi_mn_r = sqrt(eps_m*eps_n) * cos(m*pi*x/a) * cos(n*pi*y/a) ;
            psi_mn_r0 = sqrt(eps_m*eps_n) * cos(m*pi*x_0/a) * cos(n*pi*y_0/a) ;
            sum = sum + psi_mn_r.*psi_mn_r0./k_zmn .* exp(-1j.*k_zmn.*z) ;
        end
    end 
   G(ifrq) = -1j/S .* sum;
end

figure(3) ;
plot(f,20*log10((abs(G)/(2*10^-5)))) ;
grid on;
xlabel('Frequency (Hz)') ;
ylabel('Pressure level (dB SPL)') ;
title('Pressure in function of frequency') ;