%% Code to read the zigzag coordinates and velocity values and visualise them
% for velocity profile extraction
% made by AD
%clc;
close all
clear;
%% Reading, scaling and translation
tic
filename = 'Re2e5_a0p2_z0_xyzuvwp_fine_a8.txt';
surface_slice = readmatrix('surface_slice_ordered.csv');
all_line_data = load(filename);                 %all_line_data = flipud(all_line_data);
D = 2; Re = 1.4e5; mu = D/Re; %--------------------input-----------
N = size(surface_slice, 1); n1 = 500; n2 = 500; % define the number of lines and no of divisions reqd
frac1=0.15; frac2 = 1-frac1; % inner frac1 of each line will have n1 divisions
                            % outer frac2 of each line will have n2 divisions
tot_len = 0.1*D; % the height of velocity profiles are taken till 10% of Dia
intermediate_y = frac2*0 + frac1*tot_len;
y = [linspace(0,intermediate_y,n1+1)'; linspace(intermediate_y,tot_len,n2+1)'];
y(n1+1)=[];
line = cell(N,1); Tw = zeros(N,1);vplus = Tw;
theta = Tw;
P_1 = surface_slice(1,1:3)
P_2 = surface_slice(500,1:3)
P_3 = surface_slice(1000,1:3)
% Calculate vectors lying on the plane
v1 = P_2 - P_1;
v2 = P_3 - P_1;
% Calculate the normal vector
normal_vector = cross(v1, v2);
% Normalize the normal vector
normal_vector = normal_vector / norm(normal_vector);
disp('Normal vector:');
disp(normal_vector);
for i = 1:N
    theta(i) = surface_slice(i, 5);
    start = (i-1)*(n1+n2+1)+1; ending = i*(n1+n2+1);
    line{i} = all_line_data(start:ending,:);
    line{i}(:,12) = sqrt(line{i}(:,4).^2 + line{i}(:,5).^2 + line{i}(:,6).^2);   % calculating velocity magnitude
    % Calculate the norm for each row
    %norms = sqrt(line{i}(:,1).^2 + line{i}(:,2).^2 + line{i}(:,3).^2);
    line{i}(:,13:15) = [line{i}(:,1)  line{i}(:,2) line{i}(:,3)] ./ norm([line{i}(:,1)  line{i}(:,2) line{i}(:,3)]); %radial vector
    line{i}(:, 20:22) = repmat(normal_vector, size(line{i}, 1), 1);   %normal_vector % Calculate the unit vector components for each row
    line{i}(:, 23:25) = cross(line{i}(:, 20:22), line{i}(:, 13:15), 2); 
% Calculate the norm for each row of columns 23-25
    norms = sqrt(sum(line{i}(:, 23:25).^2, 2));
% Normalize each row of columns 23-25
    line{i}(:, 23:25) = line{i}(:, 23:25) ./ norms;  %tangential_vector
% Project the velocity onto the normal vector to get the normal component
   %dotProduct = sum(line{i}(:, 4:6) .* line{i}(:, 20:22), 2); 
    line{i}(:,40) = sum(line{i}(:, 4:6) .* line{i}(:, 20:22), 2); 
    % Multiply the dot product result with the vectors in columns 20-22
    line{i}(:, 26:28) = line{i}(:,40) .* line{i}(:, 20:22);  % Project the velocity onto the normal vector to get the normal component
   % line{i}(:, 26:28) = dot(line{line{i}(:,1:3), line{i}(:, 20:22)) * line{i}(:, 20:22);  %     
       
line{i}(:, 29:31) = line{i}(:, 4:6) - line{i}(:, 26:28);  % in_plane_velocity = velocity - normal_component;% Subtract the normal component from the velocity to get the in-plane velocity

%line{i}(:, 35:37) = dot(line{i}(:, 29:31), line{i}(:, 23:25),2).* line{i}(:, 23:25);  
line{i}(:, 8) = dot(line{i}(:, 29:31), line{i}(:, 23:25),2);
  % line{i}(:,7)=sqrt(line{i}(:,6).^2 + line{i}(:,7).^2);                   % total tangential speed
    %line{i}(:,9)=line{i}(:,5)*cos(theta_surf)-line{i}(:,4)*sin(theta_surf); % normal speed
    Tw_t = mu*((line{i}(2,8)-line{i}(1,8))/(y(2)-y(1))); % along x (streamwise dirn)
    cf(i)=2*Tw_t;
    Tw_w =0;% mu*((line{i}(2,6)-line{i}(1,6))/(y(2)-y(1))); % along z
    Tw(i) = sqrt(Tw_t.^2 + Tw_w.^2);%*(Tw_t./abs(Tw_t)); % vector sum of x & z component
    if Tw(i)>0; vplus(i)=sqrt(Tw(i)); elseif Tw(i)<0; vplus(i)=-sqrt(-Tw(i)); end
    line{i}(:,10)=y*vplus(i)/mu; % y+
    line{i}(:,11)=line{i}(:,8)/vplus(i); % u+
end
 
% figure(11)
% plot(Tw)


%lawp_radiantothitha=rad2deg(theta);

%lawp=180-theta
%theta = rad2deg(theta);

%theta = theta+180
%theta = 180-theta;  % this formula is correct 

%%
%seam_side =

%to_plot=[1:10:50];  % 
%to_plot=[51:10:100];  % 
%to_plot=[81:2:90] %trip 1 
%to_plot=[91:2:100] %trip 1 
%to_plot=[101:2:114] 
%to_plot=[115:2:134]   %TRIP 2
to_plot = [563]
%to_plot = [164:2:184] %trip 3
%to_plot = [185:2:205]%trip3
%to_plot = [206:2:221] 
%to_plot = [611]
%non_seam_sie
%to_plot = []
clr = {'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE';'#0072BD';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#4DBEEE'};clr2=flip(clr);
clr= ['m';clr;'r';clr2];
inner = logspace(-2,2.7782,36);
outerx= inner; outery=log(outerx)/0.41+5;
clrid = 1;
set(groot,'defaultTextInterpreter','latex')
fig=figure(1); fig.WindowState='maximized';

subplot(1,2,1);
for i=to_plot
    if clrid<=6; str='--'; elseif clrid<=12; str='-.'; elseif clrid<=18; str = ':'; else str='-'; end
    plot(line{i}(:,8),y/D,LineStyle=str,LineWidth=1.2,Color=clr{clrid},DisplayName=num2str(theta(i),5)) % num2str(i*15-15)
    hold on; clrid=clrid+1;
end
ylim([0 0.01]);
xlabel('$\bar{u}/U_{\infty} \rightarrow$',FontSize=15); ylabel('$y/D \rightarrow$',FontSize=15);
%title('Velocity profiles at arbitrary plane ($\theta = 0^\circ$ for Re=2e5,  $\alpha$=0.0 at various $\phi$)','Interpreter','latex',FontSize=18);
leg=legend(Location='north'); leg.NumColumns=8;
grid on; pbaspect([1 1 1]); %grid minor;
hold off

clrid=1;
subplot(1,2,2);
for i=to_plot
    if clrid<=6; str='--'; elseif clrid<=12; str='-.'; elseif clrid<=18; str = ':'; else str='-'; end
    semilogx(abs(line{i}(:,10)),line{i}(:,11),LineStyle=str,LineWidth=1.2,Color=clr{clrid},DisplayName=num2str(theta(i),6))
    hold on; clrid=clrid+1;
end

semilogx(inner,inner,'-.r',LineWidth=1,DisplayName='u+ = y+')
semilogx(outerx,outery,'-.b',LineWidth=1,DisplayName='u+=(ln(y+)/0.41)+5')
xlabel('$y+ \rightarrow$',FontSize=15); ylabel('$u+ \rightarrow$',FontSize=15);
title('Velocity profiles in terms of inner variables',FontSize=18);
leg=legend(Location='northwest'); leg.NumColumns=2;
grid on; pbaspect([1 1 1]); ylim([8 24]); 
xlim([0.1 600]); 
hold off
exportgraphics(gcf,'velocity_prof_2e5a0p2_a8_at_diferent theta.png',Resolution=300)

%% Plotting Cf vs Theta
fig=figure(2); fig.WindowState='maximized';
%plot (t*180/pi,cf);
% for i = 1:size(surface_slice, 1);
%     if theta(i) > 180
%         cf(i) = -cf(i);
%     end
%  end
plot (theta,cf)
xlim([0 180]);
ylabel('$C_f \rightarrow$',FontSize=15); xlabel('$\phi \rightarrow$',FontSize=15);
title('$C_f$ at arbitrary plane ($\theta = 0^\circ$)  for Re=2e5, $\alpha$=0.0 at various $\phi$','Interpreter','latex',FontSize=18); %arbitaryplane($\theta$=30  for Re=2e5 $\alpha$=0.0 at various $\phi$','Interpreter','l
leg=legend('C_f',Location='north'); leg.NumColumns=8;
grid on; %pbaspect([1 1 1]); %grid minor;
hold off
%%
l_cf = length(cf);
cf_vs_theta = zeros(l_cf,2);
cf_vs_theta(:,1) = cf; cf_vs_theta(:,2) = theta;

writematrix(cf_vs_theta,"cf_vs_theta_2e5a0p2_a8_y0.txt");

%%
figure

plot(surface_slice(:,5), surface_slice(:,4));
