%Computes junctional tension and elongation averages and their SEM
function [TensionsMean, ElongationsMean] = GetTensionsElongationsAverage(folder_path)
N = 12;
path = [folder_path 'sim_number_0/'];
[Tensions, Elongation] = GetTensionsElongations(path);
n_times = length(Tensions(:,1));

population_ratios = zeros(n_times,N);
population_hor = zeros(n_times,N);
population_vert = zeros(n_times,N);
population_all = zeros(n_times,N);
population_elongation = zeros(n_times,N);
for i=0:11
    path = [folder_path 'sim_number_' int2str(i) '/'];
    [Tensions, Elongation] = GetTensionsElongations(path);
    
    this_hor = Tensions(:,2);
    this_vert = Tensions(:,8);
    this_all = Tensions(:,11);
    ratio = this_hor./this_vert;
    population_ratios(:,i+1) = ratio;
    population_hor(:,i+1) = this_hor;
    population_vert(:,i+1) = this_vert;
    population_all(:,i+1) = this_all;
    
    this_elongation = Elongation(:,2);
    population_elongation(:,i+1) = this_elongation;
end
n_times = length(this_hor);
TensionsMean = zeros(n_times,6);
ElongationsMean = zeros(n_times,1);
for i=1:n_times
    TensionsMean(i,1) = sum(population_hor(i,:))/N;
    TensionsMean(i,2) = std(population_hor(i,:))/sqrt(N);
    TensionsMean(i,3) = sum(population_vert(i,:))/N;
    TensionsMean(i,4) = std(population_vert(i,:))/sqrt(N);
    TensionsMean(i,5) = sum(population_ratios(i,:))/N;
    TensionsMean(i,6) = std(population_ratios(i,:))/sqrt(N);
    TensionsMean(i,8) = sum(population_all(i,:))/N;
    TensionsMean(i,9) = std(population_all(i,:))/sqrt(N);
    ElongationsMean(i,1) = sum(population_elongation(i,:))/N;
    ElongationsMean(i,2) = std(population_elongation(i,:))/sqrt(N);
end

end