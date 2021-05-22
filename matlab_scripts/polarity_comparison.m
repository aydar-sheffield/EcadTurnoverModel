%Compares nematic polarity from multiple simulated conditions
close all;
clear all;
font_size = 18;
%Polarity
% folder_path = '../Test_output/TestEcadTurnover/';
folder_path = '../../../EcadSimulations/ensemble/wt_long/sim_number_0/';
[Tensions1, Elongations1] = GetTensionsElongations(folder_path);
times = Tensions1(:,1);

folder_path = '../../../EcadSimulations/ensemble/wt_long/';
[Magnitudes, Angles] = GetPolarityAverage(folder_path);
%50 unites of simulation time equals to 60 minutes
t_h = 50;
end_obs_time = (length(times)-1)/2+1;
end_obs_time = length(times);
times = times/t_h;
dn = 20;
obs_indices = 1:dn:end_obs_time;

folder_path = '../../../EcadSimulations/ensemble/p120_long/';
[Magnitudes1, Angles1] = GetPolarityAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_long/';
[Magnitudes2, Angles2] = GetPolarityAverage(folder_path);

dash_max = 1;
dash_min = -0.1;
delta = 0.01*(dash_max-dash_min);
Release_dash = dash_min:delta:dash_max;
Release_time = 3.5*ones(length(Release_dash),1);

figure(1)
p = errorbar(times(obs_indices), Magnitudes(obs_indices,1), Magnitudes(obs_indices,2), 'r','LineWidth',2);
hold on;
p1 = errorbar(times(obs_indices), Magnitudes1(obs_indices,1), Magnitudes1(obs_indices,2), 'b','LineWidth',2);
hold on;
p2 = errorbar(times(obs_indices), Magnitudes2(obs_indices,1), Magnitudes2(obs_indices,2), 'k','LineWidth',2);
hold on;
p3 = plot(Release_time, Release_dash,'.k');
legend([p,p1,p2, p3], {'wt', 'p120','dia', 'Release'},'FontAngle','italic');
xlabel('Time (h)');
ylabel('Magnitude, a.u.');
title('Long stretch, polarity magnitude');


axis([0 16 -0.1 0.9]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
box on;
saveas(gcf,'polarity_magnitude_long','fig');
saveas(gcf,'polarity_magnitude_long','png');



dash_max = 100;
dash_min = -100;
delta = 0.01*(dash_max-dash_min);
Release_dash = dash_min:delta:dash_max;
Release_time = 3.5*ones(length(Release_dash),1);
figure(2)
% H = shadedErrorBar(times(obs_indices),-Angles(obs_indices,1)*180/pi,Angles(obs_indices,2)*2*180,'lineProps','-r');
% H.mainLine.LineWidth = 2;
p = errorbar(times(obs_indices),-Angles(obs_indices,1)*180/pi,Angles(obs_indices,2)*2*180,'r', 'LineWidth',2);
hold on;
p1 = errorbar(times(obs_indices),Angles2(obs_indices,1)*180/pi,Angles2(obs_indices,2)*2*180,'k', 'LineWidth',2);
hold on;
p3 = plot(Release_time, Release_dash,'.k');
axis([0 16 -100 100]);
xlabel('Time (h)');
ylabel('Angle, degrees');
title('Long stretch, polarity angle');
legend([p, p1, p3],{'wt','dia','Release'},'FontAngle','italic');
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
box on;
saveas(gcf,'polarity_angle_long','fig');
saveas(gcf,'polarity_angle_long','png');

%%
folder_path = '../../../EcadSimulations/ensemble/wt_short/';
[Magnitudes, Angles] = GetPolarityAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_short/';
[Magnitudes1, Angles1] = GetPolarityAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_short/';
[Magnitudes2, Angles2] = GetPolarityAverage(folder_path);

dash_max = 1;
dash_min = -0.1;
delta = 0.01*(dash_max-dash_min);
Release_dash = dash_min:delta:dash_max;
Release_time = 0.5*ones(length(Release_dash),1);

figure(3)
p = errorbar(times(obs_indices), Magnitudes(obs_indices,1), Magnitudes(obs_indices,2), 'r','LineWidth',2);
hold on;
p1 = errorbar(times(obs_indices), Magnitudes1(obs_indices,1), Magnitudes1(obs_indices,2), 'b','LineWidth',2);
hold on;
p2 = errorbar(times(obs_indices), Magnitudes2(obs_indices,1), Magnitudes2(obs_indices,2), 'k','LineWidth',2);
hold on;
p3 = plot(Release_time, Release_dash,'.k');
legend([p,p1,p2, p3], {'wt', 'p120','dia', 'Release'},'FontAngle','italic');
xlabel('Time (h)');
ylabel('Magnitude, a.u.');
title('Short stretch, polarity magnitude');
axis([0 16 -0.1 0.9]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
box on;
saveas(gcf,'polarity_magnitude_short','fig');
saveas(gcf,'polarity_magnitude_short','png');

figure(4)
dash_max = 100;
dash_min = -100;
delta = 0.01*(dash_max-dash_min);
Release_dash = dash_min:delta:dash_max;
Release_time = 0.5*ones(length(Release_dash),1);
p = errorbar(times(obs_indices),-Angles(obs_indices,1)*180/pi,Angles(obs_indices,2)*2*180,'r', 'LineWidth',2);
hold on;
p1 = errorbar(times(obs_indices),Angles2(obs_indices,1)*180/pi,Angles2(obs_indices,2)*2*180,'k', 'LineWidth',2);
hold on;
p3 = plot(Release_time, Release_dash,'.k');
axis([0 16 -100 100]);
xlabel('Time (h)');
ylabel('Angle, degrees');
title('Short stretch, polarity angle');
legend([p, p2, p3],{'wt','dia','Release'},'FontAngle','italic');
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
box on;
saveas(gcf,'polarity_angle_short','fig');
saveas(gcf,'polarity_angle_short','png');

