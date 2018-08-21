function currentPerimeter = perimeter_cal(sequential_points)
% calculate the perimeter of a polygon, which is sequentially positioned by
% sequential_points
% 11-April-2016

if size(sequential_points, 1)<3
    fprintf('ERROR1 insufficient point number !\n');
end

dis = sequential_points([end 1:(end-1)], :)-sequential_points(1:(end), :);

currentPerimeter = sum(sqrt(sum(dis.^2, 2)));

end

