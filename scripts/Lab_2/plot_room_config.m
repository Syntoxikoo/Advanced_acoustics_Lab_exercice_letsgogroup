function plot_room_config(L,mics,srcs,door, filename)
    arguments
        L =[3.14, 4.38, 3.27]
        mics =  [L(1)/2, L(2)/2, L(3)/2]
        srcs = [0.16,0.155,0.155]
        door = true
        filename = "DTU_room"
    end

fig = figure;
ax = axes('Parent', fig);
view(ax, 150,20);  

daspect(ax, [1 1 1]);  
xlim(ax, [0, L(1)]);
ylim(ax, [0, L(2)]);
zlim(ax, [0, L(3)]);

% Set labels
xlabel(ax, 'x (m)');
ylabel(ax, 'y (m)');
zlabel(ax, 'z (m)');

hold on

% Plot microphones
scatter3(ax, mics(1,1), mics(1,2), mics(1,3),100, "Marker","square",'MarkerFaceColor', "black", "MarkerEdgeColor", "black" );
plot3(ax, [mics(1,1), mics(1,1)], [mics(1,2), mics(1,2)], [0, mics(1,3)], 'k-', 'LineWidth',2);
% scatter3(ax, mics(2,1), mics(2,2), mics(2,3), 'kx');

% Plot sources
scatter3(ax, srcs(1,1), srcs(1,2), srcs(1,3),100, 'r', 'o','MarkerFaceColor', "red");
plot3(ax, [srcs(1,1), srcs(1,1)], [srcs(1,2), srcs(1,2)], [0, srcs(1,3)], 'r-', 'LineWidth',2);

% Add text labels
text(ax, mics(1,1), mics(1,2), mics(1,3)+0.2, 'Mic', 'Color', 'k', 'FontSize', 12);
% text(ax, mics(2,1)-1, mics(2,2), mics(2,3), 'Mic 2-3', 'Color', 'k', 'FontSize', 10);

text(ax, srcs(1,1), srcs(1,2), srcs(1,3)+0.2, 'Src', 'Color', 'r', 'FontSize', 12);


% Define room walls
Lx = L(1);
Ly = L(2);
Lz = L(3);

vertices = [0 0 0; Lx 0 0; Lx Ly 0; 0 Ly 0];
faces = [1 2 3 4];
patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'k');

vertices = [0 0 Lz; Lx 0 Lz; Lx Ly Lz; 0 Ly Lz];
faces = [1 2 3 4];
patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'k');

vertices = [0 0 0; 0 Ly 0; 0 Ly Lz; 0 0 Lz];
faces = [1 2 3 4];
patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'k');

vertices = [Lx 0 0; Lx Ly 0; Lx Ly Lz; Lx 0 Lz];
faces = [1 2 3 4];
patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'k');

vertices = [0 Ly 0; Lx Ly 0; Lx Ly Lz; 0 Ly Lz];
faces = [1 2 3 4];
patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'k');

vertices = [0 0 0; Lx 0 0; Lx 0 Lz; 0 0 Lz];
faces = [1 2 3 4];
patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
if door == true
    vertices = [1. 0.01 0; 2 0.01 0; 2 0.01 2; 1 0.01 2];
    faces = [1 2 3 4];
    patch(ax, 'Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue', 'FaceAlpha', 0.7, 'EdgeColor', 'k');
end  

hold off; grid on

saveas(gcf, "figures/"+filename + ".eps", "epsc")