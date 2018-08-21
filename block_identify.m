function [cate_value_x, cate_value_y] = block_identify( kdtree, y_coor, x_coor, fore_x_coor, back_x_coor, fore_y_coor, back_y_coor  )
% Based on the current partitioning structure, we manually determine the
% block number for easy index



pIndex = kdtree.pIndex;
points = kdtree.points;
kdtree = kdtree.kdtree;

leafIndex = find((kdtree(:, 1)>0)&(kdtree(:, 2)==0));
% These are the added points

%% identify these blocks from Axis x
y_midy = ([0;y_coor]+[y_coor;1])/2;
y_midx = zeros(numel(y_midy), (size(fore_x_coor, 1)+1));
y_nums = zeros(1, numel(y_midy));
for ii = 1:numel(y_midy)
    forward_sel = ((y_midy(ii)-fore_y_coor).*(y_midy(ii)-back_y_coor)<0);     %%% this is for the x coordinate's calculation
    y_nums(ii) = sum(forward_sel);
    if y_nums(ii)>0
        % the cut point on the vertical direction
        scoor = fore_x_coor(forward_sel) + ((y_midy(ii)-fore_y_coor(forward_sel)).*(back_x_coor(forward_sel)-fore_x_coor(forward_sel)))./(back_y_coor(forward_sel)-fore_y_coor(forward_sel));
        y_midx(ii, 1:(y_nums(ii)+1)) = ([0 sort(scoor')]+[sort(scoor') 1])/2;
    else
        y_midx(ii, 1) = 0.5;
    end
end


% then we need to judge its position, inside which polygon
cate_value_x = zeros(size(y_midx));
for ii = 1:numel(y_nums)
    for jj = 1:(y_nums(ii)+1)
        for kk = 1:numel(leafIndex)
            kk_origin = bsxfun(@minus, points(pIndex{leafIndex(kk)}, :), [y_midx(ii, jj), y_midy(ii)]);
            if all((kk_origin([2:end 1], 1).*kk_origin(:, 2)-kk_origin(:, 1).*kk_origin([2:end 1], 2))<0)   % to check if it is empty
                cate_value_x(ii, jj) = kk;    % cate_value is the class label of each distribution
                break;
            end
        end
    end
end


%% identify these blocks from Axis y
x_midx = ([0;x_coor]+[x_coor;1])/2;
x_midy = zeros(numel(x_midx), (size(fore_y_coor, 1)+1));
x_nums = zeros(1, numel(x_midx));
for ii = 1:numel(x_midx)
    backward_sel = (((x_midx(ii)-fore_x_coor).*(x_midx(ii)-back_x_coor))<0);        %%% this is for the y coordinate's calculation
    x_nums(ii) = sum(backward_sel);
    if x_nums(ii)>0
        % the cut point on the horizonal direction
        scoor = fore_y_coor(backward_sel) + ((x_midx(ii)-fore_x_coor(backward_sel)).*(back_y_coor(backward_sel)-fore_y_coor(backward_sel)))./(back_x_coor(backward_sel)-fore_x_coor(backward_sel));
        x_midy(ii, 1:(x_nums(ii)+1)) = ([0 sort(scoor')]+[sort(scoor') 1])/2;
    else
        x_midy(ii, 1) = 0.5;
    end
end

% then we need to judge its position, inside which polygon
cate_value_y = zeros(size(x_midy));
for ii = 1:numel(x_nums)
    for jj = 1:(x_nums(ii)+1)
        for kk = 1:numel(leafIndex)
            kk_origin = bsxfun(@minus, points(pIndex{leafIndex(kk)}, :), [x_midx(ii), x_midy(ii, jj)]);
            if all((kk_origin([2:end 1], 1).*kk_origin(:, 2)-kk_origin(:, 1).*kk_origin([2:end 1], 2))<0)   % to check if it is empty
                cate_value_y(ii, jj) = kk;    % cate_value is the class label of each distribution
                break;
            end
        end
    end
end


end

