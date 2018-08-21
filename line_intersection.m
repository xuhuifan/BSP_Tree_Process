function intersect_points = line_intersection(point_1, point_2, point_3, point_4)
% calculate the intersection points of two lines
% 11-April-2016

x12 = point_1(1) - point_2(1);
x34 = point_3(1) - point_4(1);
y12 = point_1(2) - point_2(2);
y34 = point_3(2) - point_4(2);

cvalue = x12 * y34 - y12 * x34;

% if abs(cvalue) < 0.01
%     fprintf('No intersection \n');
% else
    a = point_1(1) * point_2(2) - point_1(2) * point_2(1);
    b = point_3(1) * point_4(2) - point_3(2) * point_4(1);
    
    x = (a * x34 - b * x12) / cvalue;
    y = (a * y34 - b * y12) / cvalue;
    
% end

intersect_points = [x, y];
