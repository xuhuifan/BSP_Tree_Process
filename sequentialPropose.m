function [intersects_two, intersects_index] = sequentialPropose( sequentialPoints )
% Propose a cut based on the gerative model from the paper
% 11-April-2016

% This is for the Ostomachion Process

%
% \theta should not be uniformly sampled
% WE SHOULD USE REJECTION SAMPLING TO SAMPLE \THETA

rs_correct = 1;
sampled_time = 0;
while (rs_correct)
    
    orthogonal_theta = rand*pi;
    PositionB = [cos(orthogonal_theta), sin(orthogonal_theta)];
    sequential_dis = sum(bsxfun(@times, sequentialPoints, PositionB), 2);
    
    
    seq_max = max(sequential_dis);
    seq_min = min(sequential_dis);
    largeRatio = pi*sqrt(2);
    
    if (largeRatio*rand <= (seq_max-seq_min))
        rs_correct = 0;
    end
    sampled_time = sampled_time + 1;
end
%
%
%



cut_position = seq_min*PositionB + rand*(seq_max-seq_min)*PositionB;

cut_direction = [-PositionB(2) PositionB(1)];

% calculate the projected points of sequentialPoints on the cut_position
% and cut_direction

point_3 = cut_position;
point_4 = cut_position + cut_direction;

intersects_two = [];
intersects_index = [];
for ii = 1:size(sequentialPoints, 1)
    if ii<size(sequentialPoints, 1)
        point_1 = sequentialPoints(ii, :);
        point_2 = sequentialPoints(ii+1, :);
    else
        point_1 = sequentialPoints(ii, :);
        point_2 = sequentialPoints(1, :);
    end
    intersect_point = line_intersection(point_1, point_2, point_3, point_4);
    
    if (intersect_point-point_1)*(intersect_point-point_2)'<0
        intersects_two = [intersects_two;intersect_point];
        intersects_index = [intersects_index; ii];
    end
    
end

if numel(intersects_index)~=2
    fprintf('Wrong intersection Number. \n');
end

end

