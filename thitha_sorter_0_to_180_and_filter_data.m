clear all;
clc
close all;% Load the CSV file

data = readmatrix('slice_surface.csv');

% Assume the columns are ordered as x, y, z, ...
% Filter out the rows where the y-coordinate (2nd column) is positive
filtered_data = data(data(:, 3) >= 0, :);

% Save the filtered data back to a new CSV file
writematrix(filtered_data, 'filtered_file.csv');

data = readmatrix('filtered_file.csv');
min_value = min(data(:,4));
max_value = max(data(:,4));
% Extract X, Y, Z coordinates from the data
X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

% Define the reference vector
reference_point = [-1, 0, 0];

% Function to calculate the angle between two vectors
calculate_angle = @(ref, point) acosd(dot(ref, point) / (norm(ref) * norm(point)));

% Initialize an array to store the angles
angles = zeros(size(X));

% Calculate the angle for each point
for i = 1:length(X)
    point = [X(i), Y(i), Z(i)];
    angles(i) = calculate_angle(reference_point, point);
end
% Append the angles to the original data
data_with_angles = [data, angles];
% Sort the data_with_angles by the angle column (4th column)
sorted_data_with_angles = sortrows(data_with_angles, 5);
% Write the sorted matrix to a new CSV file
writematrix(sorted_data_with_angles, 'surface_slice_ordered.csv');
