close all;
clear all;
font_size = 18;
% folder_path = '../Test_output/TestEcadTurnover/';
folder_path = '../../../EcadSimulations/ensemble/wt_long/sim_number_0/';
[Tensions1, Elongations1] = GetTensionsElongations(folder_path);
times = Tensions1(:,1);

folder_path = '../../../EcadSimulations/ensemble/wt_long/';
Ecad = GetEcadRatioAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_long/';
Ecad1 = GetEcadRatioAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_long/';
Ecad2 = GetEcadRatioAverage(folder_path);

%50 unites of simulation time equals to 60 minutes
t_h = 50;
end_obs_time = (length(times)-1)+1;
times = times/t_h;

dn = 20;

dash_max = 1.25;
dash_min = 0.75;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 3.5*ones(length(relax_dash),1);
figure(1)
p=errorbar(times(1:dn:end_obs_time), Ecad(1:dn:end_obs_time,1),Ecad(1:dn:end_obs_time,2),'r','LineWidth',2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), Ecad1(1:dn:end_obs_time,1),Ecad1(1:dn:end_obs_time,2),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), Ecad2(1:dn:end_obs_time,1),Ecad2(1:dn:end_obs_time,2),'k','LineWidth',2.5);
hold on;
p3 = plot(relax_time,relax_dash,'.k');
xlabel('Time (h)');
ylabel('Ratio');
title('Long stretch, H/V E-cadherin');
legend([p,p1, p2,p3], {'wt','p120', 'dia','Release'},'FontAngle','italic');
axis([0 16 0.75 1.25])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'Ecad_Ratio_long','fig');
saveas(gcf,'Ecad_Ratio_long','png');


folder_path = '../../../EcadSimulations/ensemble/wt_short/';
Ecad = GetEcadRatioAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_short/';
Ecad1 = GetEcadRatioAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_short/';
Ecad2 = GetEcadRatioAverage(folder_path);


dash_max = 1.25;
dash_min = 0.75;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 0.5*ones(length(relax_dash),1);
figure(2)
p=errorbar(times(1:dn:end_obs_time), Ecad(1:dn:end_obs_time,1),Ecad(1:dn:end_obs_time,2),'r','LineWidth',2.5);
hold on;
p1=errorbar(times(1:dn:end_obs_time), Ecad1(1:dn:end_obs_time,1),Ecad1(1:dn:end_obs_time,2),'b','LineWidth',2.5);
hold on;
p2=errorbar(times(1:dn:end_obs_time), Ecad2(1:dn:end_obs_time,1),Ecad2(1:dn:end_obs_time,2),'k','LineWidth',2.5);
hold on;
p3 = plot(relax_time,relax_dash,'.k');
xlabel('Time (h)');
ylabel('Ratio');
title('Short stretch, H/V E-cadherin');
legend([p,p1, p2,p3], {'wt','p120', 'dia','Release'},'FontAngle','italic');
axis([0 16 0.75 1.25])
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'Ecad_Ratio_short','fig');
saveas(gcf,'Ecad_Ratio_short','png');

