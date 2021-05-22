%This script compares simulation results with different values of U
%and constant value of rest length (i.e. independent of tension) remodelling parameter 
close all;
clear all;
% Set font size
font_size = 18;
% Load the data 
folder_path = '../../../copied_dats/k_L=10min/sim_number_0/';
[Tensions1, Elongations1] = GetTensionsElongations(folder_path);
times = Tensions1(:,1);

folder_path = '../../../copied_dats/k_L=15min_U=5/';
[Tensions, Elongations] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../copied_dats/k_L=15min_U=10/';
[Tensions1, Elongations1] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../copied_dats/k_L=15min_U=15/';
[Tensions2, Elongations2] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../copied_dats/k_L=15min_U=30/';
[Tensions3, Elongations3] = GetTensionsElongationsAverage(folder_path);

%50 unites of simulation time equals to 60 minutes
t_h = 50; 
end_obs_time = (length(times)-1)/2+1;
times = times/t_h*60;

%Data from Duda et al., 2019. 
data_recoil_ratio = [1; 3.3; 2];
data_recoil_ratio_error = [0.75/2; 0.8/2; 0.9/2];
time_recoil_ratio = [0; 24; 120];
dn = 20;
data_t_point = 1*dn+1;
data_t_point2 = 5*dn+1;

average_tensions_vert = Tensions(:,3);
average_tensions_hor = Tensions(:,1);
error_tensions_hor = Tensions(:,2);
ratio_hor_vert = Tensions(:,5);
error_ratio = Tensions(:,6);

average_tensions_vert1 = Tensions1(:,3);
average_tensions_hor1 = Tensions1(:,1);
error_tensions_hor1 = Tensions1(:,2);
ratio_hor_vert1 = Tensions1(:,5);
error_ratio1 = Tensions1(:,6);

average_tensions_vert2 = Tensions2(:,3);
average_tensions_hor2 = Tensions2(:,1);
error_tensions_hor2 = Tensions2(:,2);
ratio_hor_vert2 = Tensions2(:,5);
error_ratio2 = Tensions2(:,6);

average_tensions_vert3 = Tensions3(:,3);
average_tensions_hor3 = Tensions3(:,1);
error_tensions_hor3 = Tensions3(:,2);
ratio_hor_vert3 = Tensions3(:,5);
error_ratio3 = Tensions3(:,6);

dash_max = 7.5;
dash_min = 0.5;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 210*ones(length(relax_dash),1);
figure(1)
p = errorbar(times(1:dn:end_obs_time),ratio_hor_vert(1:dn:end_obs_time), error_ratio(1:dn:end_obs_time), 'r', 'LineWidth', 2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), ratio_hor_vert1(1:dn:end_obs_time),error_ratio1(1:dn:end_obs_time),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), ratio_hor_vert2(1:dn:end_obs_time),error_ratio2(1:dn:end_obs_time),'k','LineWidth',2.5);
hold on;
p3=errorbar(times(1:dn:end_obs_time), ratio_hor_vert3(1:dn:end_obs_time),error_ratio3(1:dn:end_obs_time),'m','LineWidth',2.5);
hold on;
p4 = plot(relax_time,relax_dash,'.k');
hold on;
pe=errorbar(time_recoil_ratio,data_recoil_ratio, data_recoil_ratio_error,'-.','LineWidth',2.5);

xlabel('Time (min)');
ylabel('Ratio');
title('Long stretch, H/V tension');
legend([pe, p,p1, p2, p3, p4], {'\textit{Exp. data}', '$U$=\textit{5}','$U$=\textit{10}',...
    '$U$=\textit{15}','$U$=\textit{30}','\textit{Release}'}, 'Interpreter', 'latex','FontAngle','italic');
axis([0 400 0.5 6.5])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'Ratio_const_rates_elasticity','fig');
saveas(gcf,'Ratio_const_rates_elasticity','png');