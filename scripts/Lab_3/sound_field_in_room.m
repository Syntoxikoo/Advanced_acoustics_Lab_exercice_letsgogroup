room =[3.28, 4.40, 3.28]; 
TR60 = 3.5;
c0 = 343;
fs=  schroeder_rect_room(room, TR60,'c0', c0);
% fs = find_f_modes(1,1,1);

lambs = c0 / fs;
N = 10000;

x =  lambs/2 + (room(1)- lambs) .*rand(N,1);
y =  lambs/2 + (room(2)- lambs) .*rand(N,1);
z =  lambs/2 + (room(3)- lambs) .*rand(N,1);

% x = room(1).*rand(N,1);
% x = room(2).*rand(N,1);
% z = zeros(N,1);
rS = [0,0,0];
rM = [x,y,z];

G = green_func_room_lab3(rM,rS,room, 'f', fs,'max_mode',10, 'absorption', true, 'T60', TR60,'no_const',true);

figure;
scatter(x,y, 10, 20*log10(abs(G)/2e-5))
axis([0 room(1) 0 room(2)])
title("p distrib accross room at schroeder freq: f=" + round(fs))
colorbar

figure;
histogram(abs(G).^2/2, 20)
title("pressure rms")

figure;
histogram(20 * log10(abs(G)/2e-5), 20)
title("pressure rms")
