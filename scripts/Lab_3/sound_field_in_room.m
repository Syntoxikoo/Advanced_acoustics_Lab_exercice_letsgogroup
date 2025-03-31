room =[3.28, 4.40, 3.28]; 
TR60 = 3.5;
c0 = 343;
fs=  schroeder_rect_room(room, TR60,'c0', c0);
% fs = find_f_modes(1,1,1);
fs = 700;
% 
% [i,j,k,f] = find_modes();
% 
% mode = [i,j,k];

lambs = c0 / fs;
M = 100;
N = 10000;
G_m =zeros(N,M);
% for i= 1: M
x =  lambs/2 + (room(1)- lambs) .*rand(N,1);
y =  lambs/2 + (room(2)- lambs) .*rand(N,1);
z =  lambs/2 + (room(3)- lambs) .*rand(N,1);
    % disp("running simulation num:"+i+"/M")

    % x = room(1).*rand(N,1);
    % y = room(2).*rand(N,1);
    % z = room(3).*rand(N,1);
rS = [0,0,0];
rM = [x,y,z];

% G= green_func_room_lab3(rM,rS,room, 'f', fs,'max_mode',10, 'absorption', true, 'T60', TR60,'no_const',true);
load("datas/Lab_3/gren_PSD1_N_10000.mat")
% end
% G = mean(G_m,1);




% figure;
% scatter(x,y, 10, 20*log10(abs(G)/2e-5))
% axis([0 room(1) 0 room(2)])
% title("p distrib accross room at schroeder freq: f=" + round(fs))
% colorbar

% stats
X_11 = abs(G).^2 /2;
X_21= 20 * log10(abs(G/sqrt(2))/2e-5);

mean_X_11 = mean(X_11);
std_X_11 = std(X_11);
rel_std_X_11 = std_X_11 / mean_X_11 ; 

mean_X_21 = mean(X_21);
std_X_21 = std(X_21);
rel_std_X_21 = std_X_21 / mean_X_21 ; 


figure;

bins = 30;
histogram(X_11, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on 
x_X11 = linspace(0,max(X_11),100);

pd_X11 = 1 / mean_X_11 * exp(-x_X11 / mean_X_11);

plot(x_X11, pd_X11, '-', 'LineWidth', 2)
title("pressure rms")

figure;

histogram(X_21, bins,'Normalization', 'pdf', 'FaceAlpha',0.4); hold on

x_X21 = linspace(max(X_21)-40,max(X_21),100);

L0 = 20 * log10(mean(abs(G)/sqrt(2))/2e-5);

pd_X21 = log(10) /10 * exp(log(10) /10 * (x_X21 - L0) - exp(log(10)/10 * (x_X21 - L0)));

plot(x_X21, pd_X21, '-', 'LineWidth', 2)
title("pressure rms")



disp("Part 4 : mean pressure : "+ round(mean_X_11,3))
disp("Part 4 : std pressure : "+ round(std_X_11,3))
disp("Part 4 : relative std pressure : "+ round(rel_std_X_11,3))


disp("Part 4 : mean pressure Level : "+ round(mean_X_21,3))
disp("Part 4 : std pressure Level : "+ round(std_X_21,3))
disp("Part 4 : relative std pressure Level : "+ round(rel_std_X_21,3))