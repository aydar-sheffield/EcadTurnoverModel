%Extracts nematic polarity and converts the results to be from -90 to 90
%degrees with respect to x-axis
function [magnitude, error_magnitude, angle, error_angle, times] = GetPolarity(path)
path0 = [path 'results_from_time_0/Polarity.dat'];
Polarity0= importdata(path0);
path1 = [path 'results_from_time_5/Polarity.dat'];
Polarity = importdata(path1);
Polarity = [Polarity0(1:end-1,:); Polarity];

times = Polarity(:,1);
magnitude = Polarity(:,2);
error_magnitude = Polarity(:,3);
angle = Polarity(:,4);
error_angle = Polarity(:,5);

for i=1:length(angle)
   if (angle(i)<=pi&&angle(i)>pi/2)
       angle(i) = (pi-angle(i));
   elseif (angle(i)>=-pi&&angle(i)<=-pi/2)
       angle(i) = angle(i)+pi;
   elseif (angle(i)>-pi/2&&angle(i)<=0)
       angle(i) = -angle(i);
   end
end



end