% Dual Axis code to plot Cp peaks and Cl vs α
% By VS
close all;
clc;
clear;
set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
ax = gca
% X-axis points
X = [0, 0.15, 0.35, 0.7, 1.0];
Cp_peak_retrating = -[-0.5156,-0.733525,-0.8072,-1.1394,-1.3177];
Cp_peak_advancing = -[-0.5061, -0.38, -0.939,-1,-0.698];
Cl = [0.0055,0.1091,-0.1604,0.081,0.3430]*10;

% Plot Cp_peak_retrating and Cp_peak_advancing on the left y-axis
yyaxis left;
set(gca, 'YColor', 'b');
hold on;
h1 = plot(X, Cp_peak_retrating, '-bo', 'LineWidth', 2.5, 'MarkerSize', 15, 'MarkerFaceColor', 'b'); % Blue circles connected by lines with filled markers
h3 = plot(X, Cp_peak_advancing, '-bs', 'LineWidth', 2.5, 'MarkerSize', 15, 'MarkerFaceColor', 'b'); % Blue squares connected by lines with filled markers

% grid on;
% grid minor;

% Add shading for the entire graph
xLimits = xlim;
yLimits = [-2 5];

% Cyan shading for x = 0 to 0.225
fill([0 0.225 0.225 0], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% Red shading for x = 0.225 to 0.55
fill([0.225 0.56 0.56 0.225], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% Yellow shading for x = 0.55 to 1
fill([0.56 1 1 0.56], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Set y-axis limits and ticks for the left axis
ylim([0 5]);
yticks(0:1:5);
ax.YMinorTick='on';
ylabel('$-\bar{C}_{P_{peak}}$', 'FontSize', 20, 'Interpreter', 'latex', 'Rotation', 90);

% Plot Cl on the right y-axis
yyaxis right;
set(gca, 'YColor', 'r');
h5 = plot(X, Cl, '-r^', 'LineWidth', 2.5, 'MarkerSize', 15, 'MarkerFaceColor', 'r'); % Red triangles connected by lines with filled markers

% Set y-axis limits and ticks for the right axis
ylim([-2 5]);
yticks(-2:1:5);
%Add texts
x = 0.1; % X-coordinate for the text
y = 5.1; % Y-coordinate for the text
text(x, y, 'M-I', 'FontSize', 20, 'Color', 'black', 'Interpreter', 'latex');
x1 = 0.375; % X-coordinate for the text

text(x1, y, 'IM', 'FontSize', 20, 'Color', 'black', 'Interpreter', 'latex');
x2 = 0.775; % X-coordinate for the text
text(x2, y, 'M-II', 'FontSize', 20, 'Color', 'black', 'Interpreter', 'latex');

ylabel('$\bar{C}_L$', 'FontSize', 20, 'Interpreter', 'latex', 'Rotation', 90);
xlabel("$ \alpha $", 'FontSize', 20);
xticks(0:0.2:1.0);
ax.XMinorTick='on';
ax.YMinorTick='on';
% Create legend with placeholders
h1 = plot(nan, nan, '-bo', 'LineWidth', 2.5, 'MarkerSize', 15, 'MarkerFaceColor', 'b');
h2 = plot(nan, nan, '-bs', 'LineWidth', 2.5, 'MarkerSize', 15, 'MarkerFaceColor', 'b');
h3 = plot(nan, nan, '-r^', 'LineWidth', 2.5, 'MarkerSize', 15, 'MarkerFaceColor', 'r');

legend([ h1, h2, h3], {'$-\bar{C}_{P_{peak}{ret}}$', '$-\bar{C}_{P_{peak}{adv}}$', '$-\bar{C}_L$',}, ...
    'Interpreter', 'latex', 'Location', 'best', 'NumColumns', 3, 'Orientation', 'vertical', ...
    'Position', [0.5, 0.80, 0.05, 0.1], 'Box', 'off');
hold off;

% Set plot aspect ratio and size
pbaspect([1 1 1]);
set(gca, 'FontSize', 20);
set(gcf, 'position', get(0, 'Screensize'));

% Export the plot as a PNG image
exportgraphics(gcf, 'Cp_peak_Ret_Adv1.png', 'Resolution', 600);
