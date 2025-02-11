%% Code to read inner endpoint,outward normal vector & create the line geometry
% for velocity profile extraction
% made by AD
% clc;
close all
clear
%% Reading, scaling and optional translation
%filename = 'xyzuvwtheta_z0_r7to9_60k_24.csv';
%filename = 'xyzuvwtheta_z01135_r7to9_60k_24.csv';
filename = 'surface_slice_ordered.csv';
spikedataraw = load(filename); % x y z XGKUN YGKUN ZGKUN, read the readme file in this dir for details
D=2; tot_len = 0.1*D; tot_len = tot_len/1; % 10% of diameter & scaling
N=length(spikedataraw(:,1));m1=500;m2=500;M=m1+m2+1; % define the number of lines and no of divisions in them
frac1=0.15; frac2 = 1-frac1; % inner frac1 of each line will have n1 divisions
                            % outer frac2 of each line will have n2 divisions
dn = tot_len*[frac1 frac2]./[m1 m2];
spikedata = zeros(N,6); line = cell(N,1); all_lines = zeros(N*M,3);
% optional translation
trnslt = 1; %input('want the orgn ofst versn lik old days,pres''0''.If u want the new corectd orgn typ,pres''1'''); % 1 will translate the coordinates such that origin lies amid top surface
if trnslt==1
    spikedata(:,1) = spikedataraw(:,1);
    spikedata(:,2) = spikedataraw(:,2);
    spikedata(:,3) = spikedataraw(:,3);
else; spikedata=spikedataraw;
end
clear spikedataraw;
for i = 1:N
    vec = [spikedata(i,1:3)]/norm([spikedata(i,1:3)]);
    temp1 = spikedata(i,1:3) + dn(1)*vec.*(0:m1)';
    temp2 = temp1(end,:)      + dn(2)*vec.*(1:m2)';
    line{i}(:,1:3)=[temp1;temp2];
    all_lines((M*(i-1)+1):(M*i),:)=line{i};
end
writematrix(all_lines,'ball_phi=0_coarsespike_1095x1000_xyz_3.csv');
% Load the existing data from the CSV file
% Load the CSV file
data = readmatrix('ball_phi=0_coarsespike_1095x1000_xyz_3.csv');

% Extract X, Y, Z coordinates from the data
X = data(:, 1);
Y = data(:, 2);
Z = data(:, 3);

% Add the header row to the data
header = {'X', 'Y', 'Z'};
data_with_header = [header; num2cell(data)];

% Write the data with header back to the CSV file
fid = fopen('ball_phi=0_coarsespike_1095x1000_xyz_3.csv', 'w');
fprintf(fid, '%s,%s,%s\n', data_with_header{1,:});
fclose(fid);
dlmwrite('ball_phi=0_coarsespike_1095x1000_xyz_3.csv', data_with_header(2:end,:), '-append', 'delimiter', ',');

%writematrix(all_lines,'xyz_z01135_r7to9_60k_spike.csv');
%writematrix(all_lines,'xyz_z0227_r7to9_60k_spike.csv');
scatter(all_lines(:,1),all_lines(:,3),Marker=".")
hold on; %daspect([1 1 1]); hold off
%exportgraphics(gcf,'frisbee_ZeroDeg_70micronOffset_1200x801_all_spikes_B_org_corect.png',Resolution=300)
