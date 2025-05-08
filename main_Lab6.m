clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))

%%
N_points = 360;
theta = linspace(0,2 * pi, N_points);
k=  1; % random for test
a= 0.1;
p0 = prad_sph_src(0,1,theta,k,a,"amp",false,"norm_ax",true);
p1 = prad_sph_src(1,1,theta,k,a,"amp",false,"norm_ax",true);
p2 = prad_sph_src(2,1,theta,k,a,"amp",false,"norm_ax",true);

plot_patternSPthe

%%

a= 0.1;
ka = 0.1;
k = ka / a;
N_points = 500;
r = logspace(log10(a), log10(1000*a), N_points);
% r = linspace(a,1000 * a, N_points);
theta0 = 0;
m = (0:2);
p = zeros(length(r),length(m));
for ii =1:length(m)
    p(:,ii) = prad_sph_src(m(ii),r,theta0,k,a,"amp",true,"norm_ax",false);
    p(:,ii) = p(:,ii)/p(end,ii);

end

plot_Lp_for_r

%% phase angle between sound pressure and particule vel

a= 0.1;
ka = 0.1;
k = ka / a;
N_points = 500;
r = logspace(log10(a), log10(1000*a), N_points);
% r = linspace(a,1000 * a, N_points);
theta0 = 0;
m = (0:2);
p = zeros(length(r),length(m));
ur = zeros(length(r),length(m));
for ii =1:length(m)
    p(:,ii) = prad_sph_src(m(ii),r,theta0,k,a,"amp",true,"norm_ax",false);
    ur(:,ii) = urad_sph_src(m(ii),r,theta0,k,a,"amp",true);
end

X = angle(p./ur);

plot_phasediff_for_r


%% task 2.1

a= 0.5;
ka_a = [0.1 1 5 10];
N_points = 500;
theta = linspace(0, 2*pi, N_points);
r = 100 *a;
m = 0:200;

p = zeros(N_points, length(ka_a));
p0 = zeros(1,4);

for ii = 1 : length(ka_a)
    ka = ka_a(ii);
    k = ka/a;
    
    p(:,ii) = prad_sph_src(m,r,theta,k,a,"amp",true,"norm_ax",false);
    p0(ii) =  prad_sph_src(m,r,0,k,a,"amp",true,"norm_ax",false);
end

plot_pattern_pS_sph

%% task 2.2

a= 0.5;
ka_a = [0.1 1 5 10];
N_points = 500;
theta = linspace(0, 2*pi, N_points);
r = 100 *a;
m = 0:100;

p = zeros(N_points, length(ka_a));
p0 = zeros(1,4);

for ii = 1 : length(ka_a)
    ka = ka_a(ii);
    k = ka/a;
    
    p(:,ii) = prad_sph_src(m,r,theta,k,a,"amp",true,"norm_ax",false, 'tricks',false);
    p0(ii) =  pS(r,k,1e-3,a);
    
    p(:,ii) = p(:,ii)/p0(ii);
end

plot_pattern_pS_sph_norm2


%%
load("cmap.mat")
N_points = 1000;
theta = linspace(0,pi, N_points/2);

phi = linspace(0,2*pi, N_points);
[theta,phi]=meshgrid(theta,phi);

figure;
set(gcf, 'Color', 'white');
grid off
M = 3;
N = 3;
i = 0;
for n = 0:N
    for m = 0:M
        if m>=n
            i = i+1;
            Yn1 = sph_harmonic(m,n,theta,phi);
            [x, y, z] = sph2cart(theta, phi-pi/2, abs(Yn1).^2);
            surf(x, y, z,"EdgeColor","none"); 
            title("Spherical Harmonic Visualization (absolute), m: "+m+", n: "+n);
            axis equal off
            colormap(smoothColormap)
            drawnow; % Force plot update
            pause(1); % Adjust speed
            frame(i) = getframe(gcf); % Capture frames for video export
        end
    end 
end


v = VideoWriter('figures/sphericalharmoabsolute.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);
writeVideo(v, frame);
close(v);

%%
funcharm(2,2)
