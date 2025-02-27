function fmax = modal_density(lx,ly,n, showPlot)

    arguments
        lx = 0.7
        ly = 1.
        n = 3
        showPlot = false
    end
c= 343;

f = linspace(0,100000,10000);
df = f(2)-f(1);
Nxy = zeros(size(f));
for ff= 1:length(f)
    Nxy(ff) = pi*lx*ly*f(ff)^2 / (c^2);
    if ff >1
        dN = Nxy(ff) - Nxy(ff-1);
        density = dN / df;
        if density > n
            disp("modale density above 10 for f ="+f(ff))
            break
        end
    end
end
fmax = round(f(ff));
if showPlot == true
    heightScale = 0.5; % Adjust height scaling if needed
    
    % ------------------------------------------------------------------------------------------------------------
    [columnwidth, ~] = get_widths();
    height = get_height() * heightScale; 
    fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
    tiled = tiledlayout(1, 1, "TileSpacing", "tight", "Padding", "loose");
    corder = colororder;
    nexttile
    plot(f(1:ff),Nxy(1:ff),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:));
    grid on;
    ylabel('Number of modes');
    xlabel('frequency');
    xlim([0 f(ff)])

end

end
