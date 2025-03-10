room =[3.14, 4.38, 3.27];
rS = [0,0,0];
x = room(1); y = room(2); z = 0;
rM = [x,y,z];
[G,~] = green_func_room(rM,rS,room, 'absorption', true); 

x = room(1)/2; y = room(2)/2; z = room(3)/2;
rM = [x,y,z];
[G2,f] = green_func_room(rM,rS,room, 'absorption', true); 

figure;
plot(f,20*log10(abs(G)/2e-5),"LineWidth",1) ;grid on; hold on
plot(f,20*log10(abs(G2)/2e-5),"LineWidth",1)
legend("S:corner - R:corner", "S: corner - R: full center")
xlim([min(f) max(f)])

