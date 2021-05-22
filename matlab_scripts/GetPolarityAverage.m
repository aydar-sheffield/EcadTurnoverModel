%Computes average polarity angles, magnitudes, and their SEM
function [Magnitudes, Angles] = GetPolarityAverage(folder_path)
N = 12;
path = [folder_path 'sim_number_0/'];
[magnitude, error_magnitude, angle, error_angle, times] = GetPolarity(path);
n_times = length(angle);

population_magnitudes = zeros(n_times,N);
population_angles = zeros(n_times,N);

for i=0:11
    path = [folder_path 'sim_number_' int2str(i) '/'];
    [magnitude, error_magnitude, angle, error_angle, times] = GetPolarity(path);
    
    population_magnitudes(:, i+1) = magnitude;
    population_angles(:, i+1) = angle;
end
Magnitudes = zeros(n_times,2);
Angles = zeros(n_times,2);
for i=1:n_times
    Magnitudes(i,1) = sum(population_magnitudes(i,:))/N;
    Magnitudes(i,2) = std(population_magnitudes(i,:))/sqrt(N);
    Angles(i,1) = sum(population_angles(i,:))/N;
    Angles(i,2) = std(population_angles(i,:))/sqrt(N);
end

end