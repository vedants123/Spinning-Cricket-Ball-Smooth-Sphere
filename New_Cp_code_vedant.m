%% Setup for LaTeX interpreter and plot formatting
% Code by Vedant
close all
clear all
clc
set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% Load and process the first dataset: 1p4_a0p0_cp_data.csv
data1 = readmatrix("1p4_a0p0_cp_data.csv");
x1 = data1(:, 1);
y1 = data1(:, 3);
theta_deg = linspace(0, 180, 100); 
theta_rad = deg2rad(theta_deg); % Convert degrees to radians
Re = 140000;
Cp = 1 - (9/4) .* sin(theta_rad).^2; % Theory


N1 = length(x1);
theta_temp1 = zeros(N1, 1);
for i = 1:N1
    theta_temp1(i) = atan2d(y1(i), x1(i));
end
data1(:, 5) = theta_temp1;

sorted1 = sortrows(data1, 5);
abstheta1 = abs(sorted1(:, 5));
[~, I1] = min(abstheta1);

p1_1 = 1; p2_1 = I1;
p3_1 = p2_1 + 1; p4_1 = length(sorted1);

sorted1(p1_1:p2_1, 5) = sorted1(p1_1:p2_1, 5) - 180;
sorted1(p1_1:p2_1, 5) = sorted1(p1_1:p2_1, 5) * -1;
sorted1(p3_1:p4_1, 5) = 180 - sorted1(p3_1:p4_1, 5);

final1 = sortrows(sorted1, 5);
theta_final1 = final1(:, 5);
Cp_final1 = final1(:, 4);

theta_0_180_1 = theta_final1(theta_final1 <= 180);
Cp_0_180_1 = Cp_final1(theta_final1 <= 180);

theta_180_360_1 = theta_final1(theta_final1 > 180);
Cp_180_360_1 = Cp_final1(theta_final1 > 180);

theta_180_360_reversed_1 = 180 - (theta_180_360_1 - 180);

%% Load and process the second dataset: 1p4_a0p15_cp_data.csv
data2 = readmatrix("1p4_a0p15_cp_data.csv");
x2 = data2(:, 1);
y2 = data2(:, 3);

N2 = length(x2);
theta_temp2 = zeros(N2, 1);
for i = 1:N2
    theta_temp2(i) = atan2d(y2(i), x2(i));
end
data2(:, 5) = theta_temp2;

sorted2 = sortrows(data2, 5);
abstheta2 = abs(sorted2(:, 5));
[~, I2] = min(abstheta2);

p1_2 = 1; p2_2 = I2;
p3_2 = p2_2 + 1; p4_2 = length(sorted2);

sorted2(p1_2:p2_2, 5) = sorted2(p1_2:p2_2, 5) - 180;
sorted2(p1_2:p2_2, 5) = sorted2(p1_2:p2_2, 5) * -1;
sorted2(p3_2:p4_2, 5) = 180 - sorted2(p3_2:p4_2, 5);

final2 = sortrows(sorted2, 5);
theta_final2 = final2(:, 5);
Cp_final2 = final2(:, 4);

theta_0_180_2 = theta_final2(theta_final2 <= 180);
Cp_0_180_2 = Cp_final2(theta_final2 <= 180);

theta_180_360_2 = theta_final2(theta_final2 > 180);
Cp_180_360_2 = Cp_final2(theta_final2 > 180);

theta_180_360_reversed_2 = 180 - (theta_180_360_2 - 180);

%% Load and process the third dataset: 1p4_a0p35.csv
data3 = readmatrix("1p4_a0p35.csv");
x3 = data3(:, 1);
y3 = data3(:, 3);

N3 = length(x3);
theta_temp3 = zeros(N3, 1);
for i = 1:N3
    theta_temp3(i) = atan2d(y3(i), x3(i));
end
data3(:, 5) = theta_temp3;

sorted3 = sortrows(data3, 5);
abstheta3 = abs(sorted3(:, 5));
[~, I3] = min(abstheta3);

p1_3 = 1; p2_3 = I3;
p3_3 = p2_3 + 1; p4_3 = length(sorted3);

sorted3(p1_3:p2_3, 5) = sorted3(p1_3:p2_3, 5) - 180;
sorted3(p1_3:p2_3, 5) = sorted3(p1_3:p2_3, 5) * -1;
sorted3(p3_3:p4_3, 5) = 180 - sorted3(p3_3:p4_3, 5);

final3 = sortrows(sorted3, 5);
theta_final3 = final3(:, 5);
Cp_final3 = final3(:, 4);

theta_0_180_3 = theta_final3(theta_final3 <= 180);
Cp_0_180_3 = Cp_final3(theta_final3 <= 180);

theta_180_360_3 = theta_final3(theta_final3 > 180);
Cp_180_360_3 = Cp_final3(theta_final3 > 180);

theta_180_360_reversed_3 = 180 - (theta_180_360_3 - 180);
%% Load 4th data set for α= 0.7
data4 = readmatrix("slice_surfacea0p7.csv");
x3 = data4(:, 1);
y3 = data4(:, 3);

N3 = length(x3);
theta_temp4 = zeros(N3, 1);
for i = 1:N3
    theta_temp4(i) = atan2d(y3(i), x3(i));
end
data4(:, 5) = theta_temp4;

sorted4 = sortrows(data4, 5);
abstheta4 = abs(sorted4(:, 5));
[~, I3] = min(abstheta4);

p1_3 = 1; p2_3 = I3;
p3_3 = p2_3 + 1; p4_3 = length(sorted4);

sorted4(p1_3:p2_3, 5) = sorted4(p1_3:p2_3, 5) - 180;
sorted4(p1_3:p2_3, 5) = sorted4(p1_3:p2_3, 5) * -1;
sorted4(p3_3:p4_3, 5) = 180 - sorted4(p3_3:p4_3, 5);

final4 = sortrows(sorted4, 5);
theta_final4 = final4(:, 5);
Cp_final4 = final4(:, 4);

theta_0_180_4 = theta_final4(theta_final4 <= 180);
Cp_0_180_4 = Cp_final4(theta_final4 <= 180);

theta_180_360_4 = theta_final4(theta_final4 > 180);
Cp_180_360_4 = Cp_final4(theta_final4 > 180);

theta_180_360_reversed_4 = 180 - (theta_180_360_4 - 180);
%% Load 4th data set for α= 1.0
data5 = readmatrix("slice_surface1p0.csv");
x3 = data5(:, 1);
y3 = data5(:, 3);

N3 = length(x3);
theta_temp5 = zeros(N3, 1);
for i = 1:N3
    theta_temp5(i) = atan2d(y3(i), x3(i));
end
data5(:, 5) = theta_temp5;

sorted5 = sortrows(data5, 5);
abstheta5 = abs(sorted5(:, 5));
[~, I3] = min(abstheta5);

p1_3 = 1; p2_3 = I3;
p3_3 = p2_3 + 1; p4_3 = length(sorted5);

sorted5(p1_3:p2_3, 5) = sorted5(p1_3:p2_3, 5) - 180;
sorted5(p1_3:p2_3, 5) = sorted5(p1_3:p2_3, 5) * -1;
sorted5(p3_3:p4_3, 5) = 180 - sorted5(p3_3:p4_3, 5);

final5 = sortrows(sorted5, 5);
theta_final5 = final5(:, 5);
Cp_final5 = final5(:, 4);

theta_0_180_5 = theta_final5(theta_final5 <= 180);
Cp_0_180_5 = Cp_final5(theta_final5 <= 180);

theta_180_360_5 = theta_final5(theta_final5 > 180);
Cp_180_360_5 = Cp_final5(theta_final5 > 180);

theta_180_360_reversed_5 = 180 - (theta_180_360_5 - 180);


%% Plotting all datasets in one figure
figure;
% Plot for the first dataset
a=plot(theta_0_180_1, Cp_0_180_1, 'Color', 'black', 'LineWidth', 1, 'LineStyle', '-');
hold on;
b=plot(theta_180_360_reversed_1, Cp_180_360_1, 'Color', 'black', 'LineWidth', 1, 'LineStyle', '--');

% Plot for the second dataset
c=plot(theta_0_180_2, Cp_0_180_2, 'Color', 'red', 'LineWidth', 1, 'LineStyle', '-');
d=plot(theta_180_360_reversed_2, Cp_180_360_2, 'Color', 'red', 'LineWidth', 1, 'LineStyle', '--');

% Plot for the third dataset
e=plot(theta_0_180_3, Cp_0_180_3, 'Color', 'blue', 'LineWidth', 1, 'LineStyle', '-');
f=plot(theta_180_360_reversed_3, Cp_180_360_3, 'Color', 'blue', 'LineWidth', 1, 'LineStyle', '--');

% Plot for the fourth dataset
g=plot(theta_0_180_4, Cp_0_180_4, 'Color', 'green', 'LineWidth', 1, 'LineStyle', '-');
h=plot(theta_180_360_reversed_4, Cp_180_360_4, 'Color', 'green', 'LineWidth', 1, 'LineStyle', '--');

% Plot for the fifth dataset
i=plot(theta_0_180_5, Cp_0_180_5, 'Color', 'magenta', 'LineWidth', 1, 'LineStyle', '-');
j=plot(theta_180_360_reversed_5, Cp_180_360_5, 'Color', 'magenta', 'LineWidth', 1, 'LineStyle', '--');
k=plot(theta_deg, Cp,'Color', "#7E2F8E", 'LineWidth', 1, 'LineStyle', '-.');
xticks(0:60:180);
xlim([0 180]);


hold on

% Set x-axis ticks and limits
xticks(0:30:180);
xlim([0 180]);
xlabel('$\theta$ ', 'Interpreter', 'latex','FontSize',16);
ylabel('$ \overline{C}_P $', 'Interpreter', 'latex','FontSize',16);
legend('α = 0.0','','α = 0.15','','α = 0.35','','α = 0.70','','α = 1.0','','Theory')
legend('FontSize', 10);
% Add grid, legend, and display
grid ("minor")
ax = gca; % Get current axes
ax.GridLineStyle = '--'; % Set grid line style to dashed
exportgraphics(gcf, 'Cp_vs_thetanew.png', 'Resolution', 600)
hold off;
