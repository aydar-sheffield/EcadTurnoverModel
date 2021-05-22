%Computes averages of tissue area and elongations and their SEM
function [Areas, Elongations] = GetTissueAreaElongationAverage(folder_path)
N = 12;
path = [folder_path 'sim_number_0/'];
[total_area, eccentricity, ratio, ratio_simple, times] = GetTissueAreaElongation(path);
n_times = length(times);


population_areas = zeros(n_times, 2);
population_elongation = zeros(n_times, 2);

for i=0:11
   path = [folder_path 'sim_number_' int2str(i) '/'];
   [total_area, eccentricity, ratio, ratio_simple, times] = GetTissueAreaElongation(path);
   area_change = (total_area/total_area(1)-1)*100;
   population_areas(:,i+1) = area_change;
   population_elongation(:,i+1) = ratio;
end

Areas = zeros(n_times,2);
Elongations = zeros(n_times,2);

for i=1:n_times
    Areas(i,1) = sum(population_areas(i,:))/N;
    Areas(i,2) = std(population_areas(i,:))/sqrt(N);
    
    Elongations(i,1) = sum(population_elongation(i,:))/N;
    Elongations(i,2) = std(population_elongation(i,:))/sqrt(N);
end


end