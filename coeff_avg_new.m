%% Code to read, plot and analyze coefficients of force and moments from multiple runs
% made by AD
% clc;
close all
clear
%% reading into cells, ND-correction & appending to arrays
N = 24; %____________________________???____________________________________
nts = zeros(N,1); dt = nts; fortrun = cell(N,1);
t = []; t_end = 0; coeff = cell(N,6); Coeff = []; D = 2;
charlen =[0.5 D]; % [charlen_in_force.F charlen_as_per_mxyz]
area = [0.19634375 0.25*pi*(D)^2]; % [area_in_force.F planform_area_as_per_mxyz]
% actual frontal area in mxyz is 0.48927 (this was for 0 deg case)

for i = 1:N
    filename = ['fortrun' num2str(i) '.21'];
    fortrun{i}= load(filename);
    nts(i) = length(fortrun{i}(:,2));
    dt(i) = fortrun{i}(1,1);
    % time wala part
    t = [t t_end+dt(i):dt(i):t_end+dt(i)*nts(i)];
    t_end = t(end);
    for j = 1:3 % processing the 6 different columns
        coeff{i,j} = fortrun{i}(:,j+1); % extraction
        % ND-correction
        coeff{i,j}=coeff{i,j}*(area(1)/area(2));
        if j>=4; coeff{i,j}=coeff{i,j}*(charlen(1)/charlen(2)); end
        if j==1; temp=[]; end; temp = [temp coeff{i,j}]; % define empty temp for each run & append
        if j==3; Coeff = [Coeff; temp]; end % after all the columns have been processed append it to main Coeff
    end
end
t = t/(D); % Non-Dimensionalisation of Time

%% Reverse mean and RMS operations on arrays
Coeff_mean = zeros(length(t),3); Coeff_rms = Coeff_mean;
Coeff_mean(end,:) = Coeff(end,:);
for i = (length(t)-1):-1:1
    for j = 1:3
        Coeff_mean(i,j)=trapz(t(i:end),Coeff(i:end,j))/(t(end)-t(i));
        Coeff_rms(i,j)=sqrt(trapz(t(i:end),((Coeff(i:end,j)-Coeff_mean(i,j)).^2))/(t(end)-t(i)));
    end
end

%% time averaging
fully_developed_after=14*D/dt(1);% input('fully developed flow is avhieved at what time step? ');
%fully_developed_after=sum(nts(1:2));
format long
Coeff_averaged = Coeff_mean(fully_developed_after,:)
format short

%% Plotting and comparing
no_of_plots = 1;
set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

if no_of_plots==1; fig=figure(1); fig.WindowState='maximized'; end
for j = 1:3
% j-th quantity, mean, RMS with time
    str = {'x';'y';'z'};%;'mx';'my';'mz'};
    if no_of_plots==1; subplot(1,3,j); elseif no_of_plots==6; fig=figure(j); end
%     plot(t(1:nts(1)),coeff{1,j},'-b','LineWidth',1,DisplayName='')
%     hold on
    for i = 1:N
        if mod(i,2)==1; c = '-b'; else; c = '-g'; end
        plot(t(sum(nts(1:i-1))+1:sum(nts(1:i))),coeff{i,j},c,'LineWidth',1,DisplayName=['$Run$ ' num2str(i)])
        hold on
    end
    plot(t,Coeff_mean(:,j),'-r',DisplayName='$Mean \leftarrow$'); hold on
    plot(t,Coeff_rms(:,j), '-m',DisplayName= '$RMS \leftarrow$'); hold on
    scatter(t(fully_developed_after),Coeff_mean(fully_developed_after,j),50,"red","filled","d",DisplayName='Avg''ing starts'); hold on
    % Plot the averaged value as a horizontal line
    yline(Coeff_averaged(j), '--k', 'LineWidth', 1, DisplayName=['$Avg = ' num2str(Coeff_averaged(j), '%.3f') '$']); hold on
    xlabel('$t U_{\infty}/D$'); ylabel(['$C_{' str{j} '}$'])
    title(['$C_{' str{j} '}$ vs ND-time $t U_{\infty}/D$'])
    if j==3; legend(Location='best',NumColumns=2); end; grid on; grid minor; pbaspect([3 3 1])
    hold off
end
exportgraphics(gcf,'Coeff_with_time.png',Resolution=600)
