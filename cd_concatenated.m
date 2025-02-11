clc;
clear;

% Import data from the CSV file
file = importdata("concatenated_a0p0.txt");
file1=importdata("cat_a0p15.txt");
file2=importdata("cat_a0p35.txt");
file3=importdata("cat_a0p7.txt");
file4=importdata("cat_a1p0.txt");
ax = gca
% Extract the timestep and CD values
timestep = file(:, 1);
S = numel(timestep);

% Non-dimensional time calculation
Time = 0.002 * (1:S)';
uinf = 1; % m/s
D = 2; % m
Time = (Time * uinf) / D;
timestep1 = file1(:, 1);
S1 = numel(timestep1);

% Non-dimensional time calculation
Time1 = 0.002 * (1:S1)';

Time1 = (Time1 * uinf) / D;

timestep2 = file2(:, 1);
S2 = numel(timestep2);

% Non-dimensional time calculation
Time2 = 0.002 * (1:S2)';

Time2 = (Time2 * uinf) / D;

timestep3 = file3(:, 1);
S3 = numel(timestep3);

% Non-dimensional time calculation
Time3 = 0.002 * (1:S3)';
uinf = 1; % m/s
D = 2; % m
Time3 = (Time3 * uinf) / D;


timestep4 = file4(:, 1);
S4 = numel(timestep4);

% Non-dimensional time calculation
Time4 = 0.002 * (1:S4)';
uinf = 1; % m/s
D = 2; % m
Time4 = (Time4 * uinf) / D;

% Extract CD values
CZ = file(:, 2)/16;

% Plot the original CD data
%figure; % Open a new figure window
plot(Time, CZ, 'k', 'DisplayName', '$C_{D}$','LineWidth',2);
hold on;

% Calculate the backward moving averages
 %Cd_bar = arrayfun(@(i) mean(CD(i:end)), 1:S)';
 %average_Cd_bar = mean(Cd_bar);

% Plot the backward moving average with a solid red line
%plot(Time, Cd_bar, '--k', 'DisplayName', 'Backward Moving Average','LineWidth',2);

% Cz
% Extract CZ values
CZ1 = file1(:, 2)/16;

% Plot the original CZ data
%figure; % Open a new figure window
plot(Time1, CZ1, 'r', 'DisplayName', '$C_{Z}$','LineWidth',2);
hold on;
CZ2 = file2(:, 2)/16;

% Plot the original CD data
%figure; % Open a new figure window
plot(Time2, CZ2, 'b', 'DisplayName', '$C_{D}$','LineWidth',2);
hold on;
CZ3 = file3(:, 2)/16;

% Plot the original CD data
%figure; % Open a new figure window
plot(Time3, CZ3, 'g', 'DisplayName', '$C_{D}$','LineWidth',2)
% Calculate the backward moving averages
%Cz_bar = arrayfun(@(i) mean(CZ(i:end)), 1:S)';
%average_Cz_bar = mean(Cz_bar);
hold on;
% Plot the backward moving average with a solid red line
%plot(Time, Cz_bar, '--b', 'DisplayName', 'Backward Moving Average','LineWidth',2);
CZ4 = file4(:, 2)/16;

% Plot the original CD data
%figure; % Open a new figure window
plot(Time4, CZ4, 'm', 'DisplayName', '$C_{D}$','LineWidth',2)
hold off;

% Set plot limits, labels, and title
ylim([0.2, 0.7]);
xlim([0,55])
xlabel(' $(tU_{\infty}/D)$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\bar{C}_{D}$', 'Interpreter', 'latex', 'FontSize', 15);
ax.XMinorTick='on';
ax.YMinorTick='on';
%title('$C_{D}$ ,$C_{Z}$ Time History for α = 0.35', 'Interpreter', 'latex', 'FontSize', 15);

% Add grid and legend
% grid on;
legend('α = 0.0','α = 0.15','α = 0.35','α = 0.70','α = 1.0') % Create the legend
 lgd.Interpreter = 'latex'; % Set LaTeX interpreter for the legend
exportgraphics(gcf, 'Cd_a0p0toa1.png', 'Resolution', 600);
% 
% 
% % Display the average value of Cd_bar
% disp(['The average value of Cd_bar is: ', num2str(average_Cd_bar)]);
% disp(['The average value of Cd_bar is: ', num2str(average_Cz_bar)]);