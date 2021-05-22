function [total_area, eccentricity, ratio, ratio_simple, times] = tissue_area_elongation(path)
Data = read_swaps_file(path);
times = Data(:,1);
total_area = Data(:,2);
n_elements = Data(:,3);
centroid_x = Data(:,4:2:end);
centroid_y = Data(:,5:2:end);

n_time = length(times);
eccentricity = zeros(n_time,1);
ratio = eccentricity;
ratio_simple = zeros(n_time,3);
for i=1:n_time
    centroid_x_clean = centroid_x(i, 1:n_elements(i))';
    centroid_y_clean = centroid_y(i, 1:n_elements(i))';
    XY = [centroid_x_clean centroid_y_clean];
    A = EllipseDirectFit(XY);
    
    a = A(1);
    b = A(2)/2;
    c = A(3);
    d = A(4)/2;
    f = A(5)/2;
    g = A(6);
    %ellipse axis
    a_p_num = 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g);
    a_p_den = (b^2-a*c)*(sqrt((a-c)^2+4*b^2)-(a+c));
    a_p = sqrt(a_p_num/a_p_den);
    
    b_p_num = 2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g);
    b_p_den = (b^2-a*c)*(-sqrt((a-c)^2+4*b^2)-(a+c));
    b_p = sqrt(b_p_num/b_p_den);
    
    if (b_p<a_p)
        eccentricity(i) =  sqrt(1-b_p^2/(a_p^2));
        ratio(i) = a_p/b_p;
    else
        eccentricity(i) =  sqrt(1-a_p^2/(b_p^2));
        ratio(i) = b_p/a_p;
    end
    
    min_x = min(centroid_x_clean);
    max_x = max(centroid_x_clean);
    min_y = min(centroid_y_clean);
    max_y = max(centroid_y_clean);
    AP = max_y - min_y;
    PD = max_x - min_x;
    ratio_simple(i,1) = PD/AP;
    ratio_simple(i,2) = PD;
    ratio_simple(i,3) = AP;
end

end