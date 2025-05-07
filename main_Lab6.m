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
    ur = urad_sph_src(m(ii),r,theta0,k,a,"amp",true);
end

X = angle(p./ur);

plot_phasediff_for_r


%% task 2.1

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
    
    p(:,ii) = prad_sph_src(m,r,theta,k,a,"amp",true,"norm_ax",false);
    p0(ii) =  prad_sph_src(m,r,0,k,a,"amp",true,"norm_ax",false);
end

figure;
for ii = 1:4
    polarplot(theta,20*log10(abs(p(:,ii)/p0(ii)) /2e-5));hold on
end
legend('ka=0.1', 'ka=1', 'ka=5', 'ka=10');

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
    
    p(:,ii) = prad_sph_src(m,r,theta,k,a,"amp",true,"norm_ax",false);
    p0(ii) =  pS(r,k,1e-3);
end

figure;
for ii = 1:4
    polarplot(theta,20*log10(abs(p(:,ii)/p0(ii)) /2e-5));hold on
end
legend('ka=0.1', 'ka=1', 'ka=5', 'ka=10');

plot_pattern_pS_sph_norm2