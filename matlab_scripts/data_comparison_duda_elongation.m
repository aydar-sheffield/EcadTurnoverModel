%This script is used to plot cell elongation and tension vs elongation
%curves
%The script was used to generate plots in Figure 4 and 6
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
end_obs_time = length(times);
times = times/t_h;
dn = 20;

average_tensions = Tensions(:,8);
average_tensions1 = Tensions1(:,8);
average_tensions2 = Tensions2(:,8);

elongation = Elongations(:, 1);
error_elongation = Elongations(:, 2);
elongation1 = Elongations1(:, 1);
error_elongation1 = Elongations1(:, 2);
elongation2 = Elongations2(:, 1);
error_elongation2 = Elongations2(:, 2);

dash_max = 4;
dash_min = 1;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 3.5*ones(length(relax_dash),1);

figure(1)
p=errorbar(times(1:dn:end_obs_time), elongation(1:dn:end_obs_time), error_elongation(1:dn:end_obs_time), 'r','LineWidth',2);
hold on;
p1=errorbar(times(1:dn:end_obs_time), elongation1(1:dn:end_obs_time), error_elongation1(1:dn:end_obs_time),'b','LineWidth',2);
hold on;
p2=errorbar(times(1:dn:end_obs_time), elongation2(1:dn:end_obs_time),error_elongation2(1:dn:end_obs_time),'k','LineWidth',2);
hold on;
p3=plot(relax_time, relax_dash,'.k');
xlabel('Time (h)');
ylabel('Elongation, a.u.');
title('Long stretch, elongation');
legend([p,p1,p2,p3], {'wt', 'p120','dia','Release'},'FontAngle','italic');
axis([0 16 1 4]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'elongation_long','fig');
saveas(gcf,'elongation_long','png');

max_elongation = max(elongation);
max_elongation1 = max(elongation1);
max_elongation2 = max(elongation2);

max_tension = max(average_tensions);
max_tension1 = max(average_tensions1);
max_tension2 = max(average_tensions2);
dn=1
Q = elongation(1:dn:end_obs_time);
sigma_hat = average_tensions(1:dn:end_obs_time);
dt = times(1+dn)-times(1);
dt = dt*60;
[K, eta, sigma, Q,Q_dot] = find_K_eta(sigma_hat, Q, dt);

Q1 = elongation1(1:dn:end_obs_time);
sigma_hat1 = average_tensions1(1:dn:end_obs_time);
[K1, eta1, sigma1, Q1,Q_dot1] = find_K_eta(sigma_hat1, Q1, dt);

Q2 = elongation2(1:dn:end_obs_time);
sigma_hat2 = average_tensions2(1:dn:end_obs_time);
[K2, eta2, sigma2, Q2, Q_dot2] = find_K_eta(sigma_hat2, Q2, dt);

leg = [' \tau_{wt} = ' num2str(floor(eta/K)) 'min'];
leg1 = [' \tau_{p120} = ' num2str(floor(eta1/K1)) 'min'];
leg2 = [' \tau_{dia} = ' num2str(floor(eta2/K2)) 'min'];

figure(2)
p = plot(elongation(1:dn:end_obs_time)/max_elongation, average_tensions(1:dn:end_obs_time)/max_tension,'r','LineWidth',2);
hold on;
p1 = plot(elongation1(1:dn:end_obs_time)/max_elongation1, average_tensions1(1:dn:end_obs_time)/max_tension1,'b','LineWidth',2);
hold on;
p2 = plot(elongation2(1:dn:end_obs_time)/max_elongation2, average_tensions2(1:dn:end_obs_time)/max_tension2,'k','LineWidth',2);
xlabel('Elongation, a.u.');
ylabel('Average tension, a.u. ');
axis([0.3 1.1 0 1.4]);
title('Long stretch, tension vs elongation');
legend([p, p1, p2],{['wt,' leg],['p120,' leg1], ['dia,' leg2]},'Location','Northwest','FontAngle','italic');
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'stress_strain_long','fig');
saveas(gcf,'stress_strain_long','png');
% 
% t_hAPF = [16 20 22 24 28 32];
% stress_data = [0.25 1.6 1.375 1 0.75 0.375];
% elongation_data = [0.05 0.1 0.2 0.225 0.25 0.125];
% stress_data_p120 = [0.25 1.25 0.875 0.75 0.6 0.5];
% elongation_data_p120 = [0.05 0.2 0.225 0.2 0.17 0.12];
% 
% normalised_stress = (stress_data-min(stress_data))/max(stress_data);
% normalised_elongation = (elongation_data-min(elongation_data))/max(elongation_data);
% normalised_stress_p120 = (stress_data_p120-min(stress_data_p120))/max(stress_data_p120);
% normalised_elongation_p120 = (elongation_data_p120-min(elongation_data_p120))/max(elongation_data_p120);
% 
% min_sim_elongation = min(elongation);
% min_sim_tension = min(average_tensions);
% min_sim_elongation_p120 = min(elongation1);
% min_sim_tension_p120 = min(average_tensions1);
% figure(5)
% plot(normalised_elongation', normalised_stress','--r');
% hold on;
% scatter(normalised_elongation, normalised_stress);
% hold on;
% plot((elongation(1:dn:end_obs_time)-min_sim_elongation)/max_elongation, (average_tensions(1:dn:end_obs_time)-min_sim_tension)/max_tension,'r','LineWidth',2);
% hold on;
% plot(normalised_elongation_p120', normalised_stress_p120','--b');
% hold on;
% scatter(normalised_elongation_p120, normalised_stress_p120);
% hold on;
% plot((elongation1(1:dn:end_obs_time)-min_sim_elongation_p120)/max_elongation1, (average_tensions1(1:dn:end_obs_time)-min_sim_tension_p120)/max_tension1,'b','LineWidth',2);

%%
folder_path = '../../../EcadSimulations/ensemble/wt_short/';
[Tensions, Elongations] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/p120_short/';
[Tensions1, Elongations1] = GetTensionsElongationsAverage(folder_path);

folder_path = '../../../EcadSimulations/ensemble/dia_short/';
[Tensions2, Elongations2] = GetTensionsElongationsAverage(folder_path);

t_h = 50;
end_obs_time = (length(times)-1)/2+1;
end_obs_time = length(times);
dn = 20;

average_tensions = Tensions(:,8);
average_tensions1 = Tensions1(:,8);
average_tensions2 = Tensions2(:,8);

elongation = Elongations(:, 1);
error_elongation = Elongations(:, 2);
elongation1 = Elongations1(:, 1);
error_elongation1 = Elongations1(:, 2);
elongation2 = Elongations2(:, 1);
error_elongation2 = Elongations2(:, 2);

dash_max = 4;
dash_min = 1;
delta = 0.01*(dash_max-dash_min);
relax_dash = dash_min:delta:dash_max;
relax_time = 0.5*ones(length(relax_dash),1);
hold off;
figure(3)
p=errorbar(times(1:dn:end_obs_time), elongation(1:dn:end_obs_time), error_elongation(1:dn:end_obs_time), 'r','LineWidth',2);
hold on;
p1=errorbar(times(1:dn:end_obs_time), elongation1(1:dn:end_obs_time), error_elongation1(1:dn:end_obs_time),'b','LineWidth',2);
hold on;
p2=errorbar(times(1:dn:end_obs_time), elongation2(1:dn:end_obs_time),error_elongation2(1:dn:end_obs_time),'k','LineWidth',2);
hold on;
p3=plot(relax_time, relax_dash,'.k');
xlabel('Time (h)');
ylabel('Elongation, a.u.');
title('Short stretch, elongation');
legend([p,p1,p2,p3], {'wt', 'p120','dia','Release'},'FontAngle','italic');
axis([0 16 1 4]);
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'elongation_short','fig');
saveas(gcf,'elongation_short','png');

max_elongation = max(elongation);
max_elongation1 = max(elongation1);
max_elongation2 = max(elongation2);

max_tension = max(average_tensions);
max_tension1 = max(average_tensions1);
max_tension2 = max(average_tensions2);

dn = 1;
Q = elongation(1:dn:end_obs_time);
sigma_hat = average_tensions(1:dn:end_obs_time);
dt = times(1+dn)-times(1);
dt = dt*60;
[K, eta, sigma, Q,Q_dot] = find_K_eta(sigma_hat, Q, dt);

Q1 = elongation1(1:dn:end_obs_time);
sigma_hat1 = average_tensions1(1:dn:end_obs_time);
[K1, eta1, sigma1, Q1,Q_dot1] = find_K_eta(sigma_hat1, Q1, dt);

Q2 = elongation2(1:dn:end_obs_time);
sigma_hat2 = average_tensions2(1:dn:end_obs_time);
[K2, eta2, sigma2, Q2, Q_dot2] = find_K_eta(sigma_hat2, Q2, dt);

leg = [' \tau_{wt} = ' num2str(floor(eta/K)) 'min'];
leg1 = [' \tau_{p120} = ' num2str(floor(eta1/K1)) 'min'];
leg2 = [' \tau_{dia} = ' num2str(floor(eta2/K2)) 'min'];

figure(4)
p = plot(elongation(1:dn:end_obs_time)/max_elongation, average_tensions(1:dn:end_obs_time)/max_tension,'r','LineWidth',2);
hold on;
p1 = plot(elongation1(1:dn:end_obs_time)/max_elongation1, average_tensions1(1:dn:end_obs_time)/max_tension1,'b','LineWidth',2);
hold on;
p2 = plot(elongation2(1:dn:end_obs_time)/max_elongation2, average_tensions2(1:dn:end_obs_time)/max_tension2,'k','LineWidth',2);
xlabel('Elongation, a.u.');
ylabel('Average tension, a.u. ');
axis([0.3 1.1 0 1.4]);
title('Short stretch, tension vs elongation');
legend([p, p1, p2],{['wt,' leg],['p120,' leg1], ['dia,' leg2]},'Location','Northwest','FontAngle','italic');
set(gca,'FontSize',font_size)
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [1000 774 713 555]);
saveas(gcf,'stress_strain_short','fig');
saveas(gcf,'stress_strain_short','png');




