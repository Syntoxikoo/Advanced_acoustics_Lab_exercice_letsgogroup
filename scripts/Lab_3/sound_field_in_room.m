room =[3.28, 4.40, 3.28]; 
TR60 = 3.5;
c0 = 343;
fs=  schroeder_rect_room(room, TR60,'c0', c0);

lambs = c0 / fs;

x =  lambs/2 + (room(1)- lambs) .*rand(100,1);
y =  lambs/2 + (room(2)- lambs) .*rand(100,1);
z =  lambs/2 + (room(3)- lambs) .*rand(100,1);
rS = [0,0,0];
rM = [x,y,z];

G = green_func_room_lab3(rM,rS,room, 'f', fs,'max_mode',10, 'absorption', true, 'T60', TR60,'no_const',true);
