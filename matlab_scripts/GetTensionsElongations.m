%Extracts junctional tension and cell elongations from test output folders
function [Tensions, Elongations] = GetTensionsElongations(folder_path)
path0 = [folder_path 'results_from_time_0/Tension.dat'];
Tensions0= importdata(path0);
path1 = [folder_path 'results_from_time_5/Tension.dat'];
Tensions = importdata(path1);
Tensions = [Tensions0(1:end-1,:); Tensions];


path1 = [folder_path 'results_from_time_0/Elongation.dat'];
Elongations0 = importdata(path1);
path1 = [folder_path 'results_from_time_5/Elongation.dat'];
Elongations = importdata(path1);
Elongations = [Elongations0(1:end-1,:); Elongations];
end