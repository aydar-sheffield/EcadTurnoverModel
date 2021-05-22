%Extracts tissue scale metrics, such as area and eccentricity
function [total_area, eccentricity, ratio, ratio_simple, times] = GetTissueAreaElongation(folder_path)
path0 = [folder_path 'results_from_time_0/AreaAndBoundary.dat'];
path = [folder_path 'results_from_time_5/AreaAndBoundary.dat'];

[total_area0, eccentricity0, ratio0, ratio_simple0, times0] = tissue_area_elongation(path0);
[total_area, eccentricity, ratio, ratio_simple, times] = tissue_area_elongation(path);

total_area = [total_area0(1:end-1); total_area];
eccentricity = [eccentricity0(1:end-1); eccentricity];
ratio= [ratio0(1:end-1); ratio];
ratio_simple = [ratio_simple0(1:end-1,:); ratio_simple];
times = [times0(1:end-1); times];

end