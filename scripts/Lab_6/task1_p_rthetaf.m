%%

clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))

c = 343;            % Speed of sound (m/s)
rho = 1.21;         % Air density (kg/m^3
a = 0.1;            % Sphere radius (m)
ka_vals = [0.1, 1, 5, 10];
modes = [0, 1, 2];

% %% 1- Plot pressure magnitude = f( theta) for m = 0, 1, 2

% N_points = 1000;
% theta = linspace(0, 2*pi, N_points);
% cth = cos(theta); % cos(theta) for Legendre

% %we only compute the polynoms becasue we normalize the pressure (which cancels the other terms)
% %so p_norm = abs(P_m(cos(theta)))/abs(P_m(1))
% P0 = legendre(0, cth); 
% P0 = P0(1,:);
% P1 = legendre(1, cth); 
% P1 = P1(1,:);
% P2 = legendre(2, cth); 
% P2 = P2(2,:);

% % normalization by axial pressure P_i(1)
% P0_norm = abs(P0) / abs(P0(1));
% P1_norm = abs(P1) / abs(P1(1));
% P2_norm = abs(P2) / abs(P2(1));

% P0_dB = 20 * log10(P0_norm/2e-5);
% P1_dB = 20 * log10(P1_norm/2e-5);
% P2_dB = 20 * log10(P2_norm/2e-5);

% % % eps = 1e-10; % avoid error maybe?
% % P0_dB = 20 * log10(P0_norm + eps);
% % P1_dB = 20 * log10(P1_norm + eps);
% % P2_dB = 20 * log10(P2_norm + eps);

% figure;
% polarplot(theta, P0_dB, 'LineWidth', 2); 
% hold on;
% polarplot(theta, P1_dB, 'LineWidth', 2);
% polarplot(theta, P2_dB, 'LineWidth', 2);
% %rlim([-40 0]);

% legend('m = 0 (Monopole)', 'm = 1 (Dipole)', 'm = 2 (Quadrupole)', 'Location', 'southoutside');
% title('Normalized Pressure Magnitude (in dB) in function of \theta (in °)');
% grid on;
%% New code

N_points = 360;
theta = linspace(0,2 * pi, N_points);
k=  1; % random for test
a= 0.1;
p0 = prad_sph_src(0,1,theta,k,a,"amp",false,"norm_ax",true);
p1 = prad_sph_src(1,1,theta,k,a,"amp",false,"norm_ax",true);
p2 = prad_sph_src(2,1,theta,k,a,"amp",false,"norm_ax",true);

plot_patternSPthe

%% 2- Pressure on Axis (theta = 0) = f(r) for ka = 0.1 (LF)

ka = 0.1;
k = ka / a;
N_points = 500;
r = logspace(log10(a), log10(1000*a), N_points);
theta0 = 0;

figure;
hold on;

% now p = Am * h_m^(2)(kr) 
for m = modes
    h = sphankel2(m, k*r);
    Pm = legendre(m, cos(theta0));
    Pm = Pm(1);
    
    p = abs(h .* Pm);
    p_norm = p / p(end);
    
    plot(r/a, 20*log10(p_norm /2e-5+ 1e-6), 'LineWidth', 2);
end

set(gca, 'XScale', 'log');
xlabel('r / a');
ylabel('Normalized Pressure (in dB)');
legend('m=0','m=1','m=2');
title('Pressure on axis in function of distance');
grid on;

%% Phase between pressure and radial particle velocity = f(r) for ka = 0.1

figure;
hold on;

for m = modes
    h = sphankel2(m, k*r);
    
    %derivative of spherical Hankel
    dh = -(m+1)./ (k*r) .* h + sphankel2(m-1, k*r);
    
    u = -1j ./ (rho * c) .* dh;
    p = h;

    phi = angle(p ./ u);
    
    plot(r/a, rad2deg(phi), 'LineWidth', 2);
end

set(gca, 'XScale', 'log');
xlabel('r / a');
ylabel('Phase Difference (in °)');
legend('m=0','m=1','m=2');
title('Phase between pressure and particle velocity in function of r/a');
grid on;

%% 4- Far-field pressure = f(theta) at ka = 0.1, 1, 5, 10 normalized by axial

N_points = 100;
theta = linspace(0, 2*pi, N_points);
cth = cos(theta);
r_far = 100 * a;

figure;
% tiledlayout(1,2)
% p2 = zeros(length(theta),length(ka_vals));

% ii = 0;
% nexttile
for ka = ka_vals
    % ii = ii +1;
    k = ka / a;
    
    p_total = zeros(size(theta));
    
    p_front = 0
    for m = 0:30 %truncated
        h_far = sphankel2(m, k*r_far);
        Pm = legendre(m, cth);
        Pm = Pm(1,:);
        dH = dr_sphankel2(m,a);
        p_total = p_total + (2*m + 1) * h_far .* Pm * dH;
        % p_front = p_front + (2*m+1 ) * h_far *dH;
        
    end
    % p_norm = p_total/p_front;
    polarplot(theta, 20*log10(p_total /2e-5), 'LineWidth', 2);hold on ;
    
end
% nexttile 
% polarplot(theta, 20*log10(abs(p2) /2e-5), 'LineWidth', 2, "LineStyle","--"); hold on;

legend('ka=0.1', 'ka=1', 'ka=5', 'ka=10');
title('Far-field Pressure in funcion of \theta (Normalized by Axial Pressure)');
%rlim([-40 0]);

%% 5- Far-field pressure normalized by free-field point source

%re-use theta from part 4
figure;

for ka = ka_vals
    k = ka / a;
    
    p_total = zeros(size(theta));
    
    for m = 0:30 %truncated
        h_far = sphankel2(m, k*r_far);
        dh_a = -(m+1)./ (k*a) .* sphankel2(m, k*a) + sphankel2(m-1, k*a);
        Pm = legendre(m, cth);
        Pm = Pm(1,:);
        
        coeff = (2*m + 1) ./ dh_a;
        p_total = p_total + coeff .* h_far .* Pm;
    end
    
    p_free = abs(exp(-1j * k * r_far) ./ r_far); %ff of point source    
    p_norm = abs(p_total) / p_free;
    
    polarplot(theta, 20*log10(p_norm + 1e-6), 'LineWidth', 2);
    hold on;
end

legend('ka=0.1', 'ka=1', 'ka=5', 'ka=10');
title('Far-field Pressure normalized by free-field point source');
%rlim([-40 20]);
