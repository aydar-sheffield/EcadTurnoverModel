%Extracts Junctional Ecad concentrations
function ecad_values = GetEcadValues(folder_path)
path0 = [folder_path 'results_from_time_0/JunctionEcad.dat'];
ecad0= importdata(path0);
path1 = [folder_path 'results_from_time_5/JunctionEcad.dat'];
ecad1 = importdata(path1);
ecad_values = [ecad0(1:end-1,:); ecad1];
end