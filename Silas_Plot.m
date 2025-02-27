% Silas Plot
clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/

% set duct dimensions
a = 1;
b = 0.7;
duct=[a, b];

% max m and n
M = 100;
N = 100;
N_modes=[N, M];

% set frequency
f1=350;
w1=2*pi*f1;
f2=150;
w2=2*pi*f2;

c=344;
rho=1.2;
p0=20e-6;

% calculate mode cutoff
mf = 2;
nf = 2;
fmn=freqMN(a,b,mf,nf)

% source vector
x0 = 0;
y0 = 0;

if x0>a
    error('x0 needs to be smaller than a!');
end

if y0>b
    error('y0 needs to be smaller than b!');
end

r0=[x0, y0, 0];



% receiver vector
x= a;
y= b;
rec=[x, y];
z= linspace(0,10,1000);


% Calculate Green's function

Gz1 = Gf_duct(rec,r0,z,duct,f1,N_modes,[],1,false);
Gz2 = Gf_duct(rec,r0,z,duct,f2,N_modes,[],1,false);



LGz1=20*log10(abs(Gz1)/p0);
LGz2=20*log10(abs(Gz2)/p0);

figure;
plot(z, LGz1,LineWidth=2,LineStyle="-",Color="0 0.4470 0.7410");
hold on;
plot(z, LGz2,LineWidth=2,LineStyle="--",Color="0 0.4470 0.7410");
grid on;
hold off;
set(gca, 'GridLineStyle', '-', 'LineWidth', 1.5); % Thicker grid lines and axis lines
ax = gca;
ax.XAxis.LineWidth = 2; % Thicker X-axis line
ax.YAxis.LineWidth = 2; % Thicker Y-axis line
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlabel('distance in m', 'FontSize', 14); % Label for X-axis
ylabel('Sound pressure in dBSPL', 'FontSize', 14); % Label for Y-axis
title('Sound pressure over distances', 'FontSize', 16); % Big title
legend(['f = ', num2str(f1), 'Hz'],['f = ', num2str(f2), 'Hz'], FontSize=14);


% across different frequencies
% source vector
x0 = 0;
y0 = 0;

if x0>a
    error('x0 needs to be smaller than a!');
end

if y0>b
    error('y0 needs to be smaller than b!');
end

r0=[x0, y0, 0];


% receiver vector
x= a;
y= b;
rec=[x, y];
z1= 0.1;
z2= 10;
r1=[x, y, z1];
r2=[x, y, z2];

f=linspace(0,400,1000);




Gf1 = Gf_duct(rec,r0,z1,duct,f,N_modes,[],1,false);
Gf2 = Gf_duct(rec,r0,z2,duct,f,N_modes,[],1,false);



LGf1=20*log10(abs(Gf1)/p0);
LGf2=20*log10(abs(Gf2)/p0);

figure;
plot(f, LGf1,LineWidth=2,LineStyle="-",Color="0 0.4470 0.7410");
hold on;
plot(f, LGf2,LineWidth=2,LineStyle="--",Color="0 0.4470 0.7410");
hold off;
grid on;
set(gca, 'GridLineStyle', '-', 'LineWidth', 1.5); % Thicker grid lines and axis lines
ax = gca;
ax.XAxis.LineWidth = 2; % Thicker X-axis line
ax.YAxis.LineWidth = 2; % Thicker Y-axis line
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlabel('frequency in Hz', 'FontSize', 14); % Label for X-axis
ylabel('Sound pressure in dBSPL', 'FontSize', 14); % Label for Y-axis
title('Sound pressure over frequency', 'FontSize', 16); % Big title
legend(['z = ', num2str(z1), 'm'],['z = ', num2str(z2), 'm'], FontSize=14);

