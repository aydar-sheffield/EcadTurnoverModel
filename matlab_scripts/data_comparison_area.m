%This script is used to analyse tissue scale data: area and elongation
%It was used to generate plots on Figure 6
close all;
clear all;
font_size = 18;

folder_path = '../../../EcadSimulations/ensemble/wt_long/sim_number_0/';
[Tensions1, Elongations1] = GetTensionsElongations(folder_path);
times = Tensions1(:,1);

folder_path = '../../../EcadSimulations/ensemble/wt_long/';
[Areas, Elongations] = GetTissueAreaElongationAverage(folder_path);

%50 unites of simulation time equals to 60 minutes
t_h = 50;
end_obs_time = (length(times)-1)/2+1;
times = times/t_h;
dn = 20;

folder_path = '../../../EcadSimulations/ensemble/p120_long/';
[Areas1, Elongations1] = GetTissueAreaElongationAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_long/';
[Areas2, Elongations2] = GetTissueAreaElongationAverage(folder_path);

dash_max = 100;
dash_min = -15;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 3.5*ones(length(relax_dash),1);

figure(1)
p=errorbar(times(1:dn:end),Areas(1:dn:end,1), Areas(1:dn:end,2), 'r','LineWidth',2);
hold on;
p1=errorbar(times(1:dn:end),Areas1(1:dn:end,1), Areas1(1:dn:end,2),'b','LineWidth',2);
hold on;
p2=errorbar(times(1:dn:end),Areas2(1:dn:end,1), Areas2(1:dn:end,2),'k','LineWidth',2);
hold on;
p3=plot(relax_time, relax_dash,'.k');
legend([p,p1,p2, p3], {'wt', 'p120','dia', 'Release'},'FontAngle','italic');
xlabel('Time (h)');
ylabel('Area change, %');
title('Long stretch, tissue area change');
axis([0 16 -15 100]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'tissue_area_long','fig');
saveas(gcf,'tissue_area_long','png');

dash_max = 4;
dash_min = 1;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 3.5*ones(length(relax_dash),1);
figure(2)
p=errorbar(times(1:dn:end),Elongations(1:dn:end,1), Elongations(1:dn:end,2), 'r','LineWidth',2);
hold on;
p1=errorbar(times(1:dn:end),Elongations1(1:dn:end,1), Elongations1(1:dn:end,2),'b','LineWidth',2);
hold on;
p2=errorbar(times(1:dn:end),Elongations2(1:dn:end,1), Elongations2(1:dn:end,2),'k','LineWidth',2);
hold on;
p3=plot(relax_time, relax_dash,'.k');
legend([p,p1,p2, p3], {'wt', 'p120','dia', 'Release'},'FontAngle','italic');
xlabel('Time (h)');
ylabel('Ratio, major/minor');
title('Long stretch, major/minor axis ratio');
axis([0 16 1 4]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'tissue_major_minor_long','fig');
saveas(gcf,'tissue_major_minor_long','png');

%%
folder_path = '../../../EcadSimulations/ensemble/wt_short/';
[Areas, Elongations] = GetTissueAreaElongationAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_short/';
[Areas1, Elongations1] = GetTissueAreaElongationAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_short/';
[Areas2, Elongations2] = GetTissueAreaElongationAverage(folder_path);

dash_max = 100;
dash_min = -15;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 0.5*ones(length(relax_dash),1);
hold off;

figure(3)
p=errorbar(times(1:dn:end),Areas(1:dn:end,1), Areas(1:dn:end,2), 'r','LineWidth',2);
hold on;
p1=errorbar(times(1:dn:end),Areas1(1:dn:end,1), Areas1(1:dn:end,2),'b','LineWidth',2);
hold on;
p2=errorbar(times(1:dn:end),Areas2(1:dn:end,1), Areas2(1:dn:end,2),'k','LineWidth',2);
hold on;
p3=plot(relax_time, relax_dash,'.k');
legend([p,p1,p2, p3], {'wt', 'p120','dia', 'Release'},'FontAngle','italic');
xlabel('Time (h)');
ylabel('Area change, %');
title('Short stretch, tissue area change');
axis([0 16 -15 100]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'tissue_area_short','fig');
saveas(gcf,'tissue_area_short','png');

dash_max = 4;
dash_min = 1;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 0.5*ones(length(relax_dash),1);
hold off;
figure(4)
p=errorbar(times(1:dn:end),Elongations(1:dn:end,1), Elongations(1:dn:end,2), 'r','LineWidth',2);
hold on;
p1=errorbar(times(1:dn:end),Elongations1(1:dn:end,1), Elongations1(1:dn:end,2),'b','LineWidth',2);
hold on;
p2=errorbar(times(1:dn:end),Elongations2(1:dn:end,1), Elongations2(1:dn:end,2),'k','LineWidth',2);
hold on;
p3=plot(relax_time, relax_dash,'.k');
legend([p,p1,p2, p3], {'wt', 'p120','dia', 'Release'},'FontAngle','italic');
xlabel('Time (h)');
ylabel('Ratio, major/minor');
title('Short stretch, major/minor axis ratio');
axis([0 16 1 4]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'tissue_major_minor_short','fig');
saveas(gcf,'tissue_major_minor_short','png');

