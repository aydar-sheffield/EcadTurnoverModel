%This script is used to compare junctional tension, and its anisotropy
%from short and long stretches and in different genetic backgrounds
%This script was used to generate results in Figure 4
close all;
clear all;
font_size = 18;

folder_path = '../../../EcadSimulations/ensemble/wt_long/sim_number_0/';
[Tensions1, Elongations1] = GetTensionsElongations(folder_path);
times = Tensions1(:,1);

folder_path = '../../../EcadSimulations/ensemble/wt_long/';
[Tensions, Elongations] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_long/';
[Tensions1, Elongations1] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_long/';
[Tensions2, Elongations2] = GetTensionsElongationsAverage(folder_path);

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

dash_max = 6.5;
dash_min = 0.5;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 210*ones(length(relax_dash),1);
figure(1)
pe=errorbar(time_recoil_ratio,data_recoil_ratio, data_recoil_ratio_error,'-.','LineWidth',2.5);
hold on;
% H = shadedErrorBar(times(1:dn:end_obs_time),ratio_hor_vert(1:dn:end_obs_time),error_ratio(1:dn:end_obs_time),'lineProps','-r');
% p = H.mainLine;
% p.LineWidth = 2;
% p = H.mainLine;
p = errorbar(times(1:dn:end_obs_time),ratio_hor_vert(1:dn:end_obs_time), error_ratio(1:dn:end_obs_time), 'r', 'LineWidth', 2.5);
% p=plot(times(1:dn:end_obs_time), ratio_hor_vert(1:dn:end_obs_time),'r','LineWidth',2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), ratio_hor_vert1(1:dn:end_obs_time),error_ratio1(1:dn:end_obs_time),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), ratio_hor_vert2(1:dn:end_obs_time),error_ratio2(1:dn:end_obs_time),'k','LineWidth',2.5);
hold on;
p3 = plot(relax_time,relax_dash,'.k');
xlabel('Time (min)');
ylabel('Ratio');
title('Long stretch, H/V tension');
legend([pe, p,p1, p2,p3], {'Exp. data', 'wt','p120', 'dia','Release'},'FontAngle','italic');
axis([0 400 0.5 6.5])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'Ratio_long','fig');
saveas(gcf,'Ratio_long','png');

relax_dash = 0:0.01*3.5:3.5;
relax_time = 210*ones(length(relax_dash),1);
figure(3)
p=errorbar(times(1:dn:end_obs_time), average_tensions_hor(1:dn:end_obs_time), error_tensions_hor(1:dn:end_obs_time), 'r','LineWidth',2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), average_tensions_hor1(1:dn:end_obs_time), error_tensions_hor1(1:dn:end_obs_time),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), average_tensions_hor2(1:dn:end_obs_time), error_tensions_hor2(1:dn:end_obs_time),'k','LineWidth',2.5);
hold on;
p3 = plot(relax_time,relax_dash,'.k');
xlabel('Time (min)');
ylabel('Tension, a.u.');
title('Long stretch, horizontal tension');
legend([p,p1, p2,p3], {'wt','p120', 'dia','Release'},'FontAngle','italic');
axis([0 400 0 3.5])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'horizontal_long','fig');
saveas(gcf,'horizontal_long','png');


folder_path = '../../../EcadSimulations/ensemble/wt_short/';
[Tensions, Elongations] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_short/';
[Tensions1, Elongations1] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_short/';
[Tensions2, Elongations2] = GetTensionsElongationsAverage(folder_path);

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

dash_max = 6.5;
dash_min = 0.5;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 30*ones(length(relax_dash),1);
figure(2)
p = errorbar(times(1:dn:end_obs_time),ratio_hor_vert(1:dn:end_obs_time), error_ratio(1:dn:end_obs_time), 'r', 'LineWidth', 2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), ratio_hor_vert1(1:dn:end_obs_time),error_ratio1(1:dn:end_obs_time),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), ratio_hor_vert2(1:dn:end_obs_time),error_ratio2(1:dn:end_obs_time),'k','LineWidth',2.5);
hold on;
p3 = plot(relax_time,relax_dash,'.k');
xlabel('Time (min)');
ylabel('Ratio');
title('Short stretch, H/V tension');
legend([p,p1, p2,p3], {'wt','p120', 'dia','Release'},'FontAngle','italic');
axis([0 400 0.5 6.5])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'Ratio_short','fig');
saveas(gcf,'Ratio_short','png');


relax_dash = 0:0.01*3.5:3.5;
relax_time = 30*ones(length(relax_dash),1);
figure(4)
p=errorbar(times(1:dn:end_obs_time), average_tensions_hor(1:dn:end_obs_time), error_tensions_hor(1:dn:end_obs_time), 'r','LineWidth',2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), average_tensions_hor1(1:dn:end_obs_time), error_tensions_hor1(1:dn:end_obs_time),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), average_tensions_hor2(1:dn:end_obs_time), error_tensions_hor2(1:dn:end_obs_time),'k','LineWidth',2.5);
hold on;
p3 = plot(relax_time,relax_dash,'.k');
xlabel('Time (min)');
ylabel('Tension, a.u.');
title('Short stretch, horizontal tension');
legend([p,p1, p2,p3], {'wt','p120', 'dia','Release'},'FontAngle','italic');
axis([0 400 0 3.5])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'horizontal_short','fig');
saveas(gcf,'horizontal_short','png');
