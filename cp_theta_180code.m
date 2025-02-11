clc;
clear all;
close all;

% List of input files
fileList = {'surface_slice_ordered0p0.csv', 'surface_slice_ordered0p15.csv', ...
            'surface_slice_ordered0p35.csv', 'surface_slice_ordered0p7.csv', ...
            'surface_slice_ordered1p0.csv'};

% Define colors for the plots
colors = {'k', 'r', 'b', 'g', 'm'};

% Create a figure
figure;
hold on;

% Loop through each file and plot the results
for i = 1:length(fileList)
    % Load data from the CSV file
    data = readmatrix(fileList{i});

    % Extract the relevant columns (assuming Cp is the 4th column)
    Cp = data(:, 4);

    % Define theta from 0 to 180 degrees
    theta = linspace(0, 180, numel(Cp));

    % Compute Cp(180 - theta)
    Cp_180_theta = flip(Cp);

    % Compute the desired function: [Cp(theta) - Cp(180 - theta)] * cos(theta)
    result = (Cp - Cp_180_theta) .* cosd(theta);

    % Extract the portion for plotting (0 to 90 degrees)
    theta_plot = theta(theta <= 90);
    result_plot = result(1:length(theta_plot));

    % Plot the result
    plot(theta_plot, result_plot, 'LineWidth', 2, 'Color', colors{i}, ...
        'DisplayName', sprintf('File %d', i));
end

% Label the plot
xlabel('$\theta (degrees)$');
ylabel('$[\bar{C}_p(\theta) - \bar{C}_p(180 - \theta)] cos(\theta)$');
%title('Plot of [C_p(\theta) - C_p(180 - \theta)] cos(\theta) vs \theta');
legend('α = 0','α = 0.15','α = 0.35','α = 0.70','α = 1.0')
grid on;
exportgraphics(gcf, 'cp_180-theta.png', 'Resolution', 600);


% Add a legend
legend show;
hold off;
