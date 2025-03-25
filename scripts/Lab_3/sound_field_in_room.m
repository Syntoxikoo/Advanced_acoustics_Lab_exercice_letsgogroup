room =[3.28, 4.40, 3.28]; 
TR60 = 3.5;
c0 = 343;
fs=  schroeder_rect_room(room, TR60,'c0', c0);
% fs = find_f_modes(1,1,1);
% fs = 700;

lambs = c0 / fs;
N = 1000;

x =  lambs/2 + (room(1)- lambs) .*rand(N,1);
y =  lambs/2 + (room(2)- lambs) .*rand(N,1);
z =  lambs/2 + (room(3)- lambs) .*rand(N,1);

% x = room(1).*rand(N,1);
% y = room(2).*rand(N,1);
% z = room(3).*rand(N,1);
rS = [0,0,0];
rM = [x,y,z];

G = green_func_room_lab3(rM,rS,room, 'f', fs,'max_mode',10, 'absorption', true, 'T60', TR60,'no_const',true);

figure;
scatter(x,y, 10, 20*log10(abs(G)/2e-5))
axis([0 room(1) 0 room(2)])
title("p distrib accross room at schroeder freq: f=" + round(fs))
colorbar

figure;
X_11 = abs(G);
bins = 30;
histogram(X_11, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on 
x = linspace(0,max(X_11),100);
% pd = fitdist(X_11, 'Beta');
% y = pdf(pd, x);
% plot(x, y, 'r-', 'LineWidth', 2)
title("pressure rms")

figure;
X_21= 20 * log10(abs(G)/2e-5);
histogram(X_21, bins,'Normalization', 'pdf', 'FaceAlpha',0.4); hold on

x = linspace(0,max(X_21),100);
pd = fitdist(X_21, 'Normal');
y = pdf(pd, x);
plot(x, y, 'r-', 'LineWidth', 2)
title("pressure rms")


% stats
mean_X_11 = mean(X_11);
std_X_11 = std(X_11);
rel_std_X_11 = std_X_11 / mean_X_11 * 100; 

mean_X_21 = mean(X_21);
std_X_21 = std(X_21);
rel_std_X_21 = std_X_21 / mean_X_21 * 100; 

disp("Part 4 : mean pressure : "+ round(mean_X_11,3))
disp("Part 4 : std pressure : "+ round(std_X_11,3))
disp("Part 4 : relative std pressure : "+ round(rel_std_X_11,3))


disp("Part 4 : mean pressure Level : "+ round(mean_X_21,3))
disp("Part 4 : std pressure Level : "+ round(std_X_21,3))
disp("Part 4 : relative std pressure Level : "+ round(rel_std_X_21,3))