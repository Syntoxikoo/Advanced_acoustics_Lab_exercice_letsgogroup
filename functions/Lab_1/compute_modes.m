function fmn = compute_modes(lx,ly,N_modes)
arguments
    lx double
    ly double
    N_modes =[10,10];
end
[nx, ny] = ndgrid((0:N_modes(1)),(0:N_modes(2)));
fmn = 343 / 2 *sqrt((nx/lx).^2 + (ny/ly).^2);
end