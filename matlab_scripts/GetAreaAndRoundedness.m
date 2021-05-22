%Computes tissue area and eccentricity
function [total_areas, wing_areas, eccentricity,ratio] = GetAreaAndRoundedness(path)
Data = read_swaps_file(path);
times = Data(:,1);
total_areas = Data(:,2);
wing_areas = Data(:,3);
hinge_areas = Data(:,4);
n_elements = Data(:,5);
centroid_x = Data(:,6:2:end);
centroid_y = Data(:,7:2:end);

n_time = length(times);
centroid_x_clean = centroid_x(n_time, 1:n_elements(n_time))';
centroid_y_clean = centroid_y(n_time, 1:n_elements(n_time))';
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

eccentricity = 0;
ratio = 0;
if (b_p<a_p)
    eccentricity =  sqrt(1-b_p^2/(a_p^2));
    ratio = a_p/b_p;
else
    eccentricity =  sqrt(1-a_p^2/(b_p^2));
    ratio = b_p/a_p;
end

%Convert the A to str
a = num2str(A(1));
b = num2str(A(2));
c = num2str(A(3));
d = num2str(A(4));
e = num2str(A(5));
f = num2str(A(6));

%Equation
eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
xmin=min(centroid_x_clean);
xmax=1.2*max(centroid_x_clean);
if (xmin<0)
    xmin = 1.4*xmin;
end
y_min = min(centroid_y_clean);
y_max = max(centroid_y_clean);

figure
k = boundary(centroid_x_clean, centroid_y_clean);
plot(centroid_x_clean, centroid_y_clean,'.');
hold on;
plot(centroid_x_clean(k), centroid_y_clean(k));
hold on;
ez = ezplot(eqt,[xmin xmax]);
axis([xmin xmax 1.1*min(y_min,min(ez.YData)) 1.1*max(y_max,max(ez.YData))]);
title(['']);
pbaspect([1 1 1])


n_time = 1;
centroid_x_clean = centroid_x(n_time, 1:n_elements(n_time))';
centroid_y_clean = centroid_y(n_time, 1:n_elements(n_time))';
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

eccentricity = 0;
ratio = 0;
if (b_p<a_p)
    eccentricity =  sqrt(1-b_p^2/(a_p^2));
    ratio = a_p/b_p;
else
    eccentricity =  sqrt(1-a_p^2/(b_p^2));
    ratio = b_p/a_p;
end

%Convert the A to str
a = num2str(A(1));
b = num2str(A(2));
c = num2str(A(3));
d = num2str(A(4));
e = num2str(A(5));
f = num2str(A(6));

%Equation
eqt= ['(',a, ')*x^2 + (',b,')*x*y + (',c,')*y^2 + (',d,')*x+ (',e,')*y + (',f,')'];
xmin=min(centroid_x_clean);
xmax=1.2*max(centroid_x_clean);
if (xmin<0)
    xmin = 1.4*xmin;
end
y_min = min(centroid_y_clean);
y_max = max(centroid_y_clean);

hold on;
k = boundary(centroid_x_clean, centroid_y_clean);
plot(centroid_x_clean, centroid_y_clean,'.');
hold on;
plot(centroid_x_clean(k), centroid_y_clean(k));
hold on;
ez = ezplot(eqt,[xmin xmax]);
% axis([xmin xmax 1.1*min(y_min,min(ez.YData)) 1.1*max(y_max,max(ez.YData))]);
title(['']);
pbaspect([1 1 1])
end