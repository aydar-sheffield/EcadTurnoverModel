%Extract values of Ecad and plots of SEM
function Ecad = GetEcadRatioAverage(folder_path)
N = 12;
path = [folder_path 'sim_number_0/'];
[Tensions, Elongation] = GetTensionsElongations(path);
n_times = length(Tensions(:,1));

population_ratio = zeros(n_times,N);

for i=0:11
    path = [folder_path 'sim_number_' int2str(i) '/'];
    Ecads = GetEcadValues(path);
    this_vert = Ecads(:,8);
    this_hor = Ecads(:,2);
    this_ratio = this_hor./this_vert;
    population_ratio(:,i+1) = this_ratio;
end
Ecad = zeros(n_times,2);
for i=1:n_times
    Ecad(i,1) = sum(population_ratio(i,:))/N;
    Ecad(i,2) = std(population_ratio(i,:))/sqrt(N);
end

end