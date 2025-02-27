function fmax = modal_density(lx,ly,n, plot)

    arguments
        lx = 0.7
        ly = 1.
        n = 10
        plot = false
    end
c= 343;

f = linspace(0,10000,1000);
Nxy = zeros(size(f));
for ff= 1:length(f)
    Nxy(ff) = pi*lx*ly*f(ff)^2 / (c^2);
    if ff >1
        density = Nxy(ff)-Nxy(ff-1) / (f(2)-f(1));
        if density > n
            disp("modale density above 10 for f ="+f(ff))
            break
        end
    end
end
fmax = round(f(ff));
if plot == true
    figure;
    plot(Nxy)
end

end
