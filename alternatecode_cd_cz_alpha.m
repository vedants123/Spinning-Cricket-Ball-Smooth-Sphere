clear all
close all
clc

% Data
alpha = [0.0, 0.15, 0.35, 0.7, 1.0]; % Replace with your actual alpha values
Cz = [0.0055,0.1091,-0.1604,0.081,0.3430];   % Replace with your actual Cz values
Cd = [0.4643,0.470,0.3639,0.252,0.3608]; % Replace with your actual Cd values
markerColors = {'black', 'red', 'blue', 'green', 'magenta'}; % Marker colors
set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

% Plotting
figure;
hold on;

% Define shaded regions
fill([0 0.225 0.225 0], [-1 -1 1 1], 'cyan', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([0.225 0.573 0.573 0.225], [-1 -1 1 1], 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([0.573 1 1 0.573], [-1 -1 1 1], 'yellow', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot Cz with dashed line and colored markers
for i = 1:length(alpha)
    plot(alpha(i), Cz(i), 'o', 'LineWidth', 1.5, 'MarkerSize', 10, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', markerColors{i});
end
plot(alpha, Cz, '-.', 'LineWidth', 1.5, 'Color', 'black', 'DisplayName', '${\overlineC}_Z$');

% Plot Cd with solid line and colored markers
for i = 1:length(alpha)
    plot(alpha(i), Cd(i), 's', 'LineWidth', 1.5, 'MarkerSize', 10, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', markerColors{i});
end
plot(alpha, Cd, '-', 'LineWidth', 1.5, 'Color', 'black', 'DisplayName', '${\overlineC}_D$');

% Add horizontal line at y=0
yline(0, ':k', 'LineWidth', 1);

% Get the axes object and set minor ticks
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

% Add text annotations
xText1 = 0.110; % Horizontal position (center of the plot)
yText1 = 0.51; % Vertical position (near the top)
text(xText1, yText1, 'M-I', 'HorizontalAlignment', 'center', 'FontSize', 10);
xText2 = 0.40; % Horizontal position (center of the plot)
yText2 = 0.51; % Vertical position (near the top)
text(xText2, yText2, 'IM', 'HorizontalAlignment', 'center', 'FontSize', 10);
xText3 = 0.80; % Horizontal position (center of the plot)
yText3 = 0.51; % Vertical position (near the top)
text(xText3, yText3, 'M-II', 'HorizontalAlignment', 'center', 'FontSize', 10);

% Add labels
xlabel('$$\alpha$$', 'FontSize', 10, 'Interpreter', 'latex');
ylabel('Coefficients', 'FontSize', 10);

% Add legend
legend('','','','','','','','','$\overline{C}_Z$','','','','','', '$\overline{C}_D$', 'Interpreter', 'latex', 'Location', 'northeast');

% Set axis limits for better visualization
xticks(0:0.2:1);
ylim([-0.2, 0.5]);
pbaspect([1 1 1]);

hold off;

% Save the figure with high resolution and 1:1 aspect ratio
exportgraphics(gcf, 'HighQualityPlot_1to1new.png', 'Resolution', 600);

