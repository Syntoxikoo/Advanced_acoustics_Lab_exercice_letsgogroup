room =[3.14, 4.38, 3.27];
rS = [0,0,0];
x = linspace(0,room(1),10);
y = linspace(0,room(2),10);
z = ones(size(x)) * room(3);
rM = [x',y',z'];
[G,f] = green_func_room(rM,rS,room, 'absorption', true);

idx = find_f_modes(1,1,1);
contourf(x,y,abs(G(:,:,idx)))
title("this is an example of the green function from simulation")