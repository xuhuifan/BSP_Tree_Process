function z_label_seq = z_label_determine(currentParticle, coor_eta, coor_xi, points, kdtree, nodeNum, dataNum)


%
%
% model paramters

z_label_seq = cell((nodeNum-1)/2, 1);
z_label_seq{1} = ones(dataNum, dataNum);

%% First, we need to check out the available points
for dd =1:length(z_label_seq)
    addIndex = kdtree(2:nodeNum, 1)+3;   % addIndex denote the index of the TRUE added points

    addPoints = points(addIndex, :);  % addedPoints denote the coordiantes of the added points

    fore_x_coor = addPoints((1:2:end), 1);
    back_x_coor = addPoints((2:2:end), 1);
    fore_y_coor = addPoints((1:2:end), 2);
    back_y_coor = addPoints((2:2:end), 2);

    y_coor = addPoints(:, 2);
    x_coor = addPoints(:, 1);

    y_coor(y_coor == 0) = [];
    y_coor(y_coor == 1) = [];
    x_coor(x_coor == 0) = [];
    x_coor(x_coor == 1) = [];
    y_coor = sort(y_coor);
    x_coor = sort(x_coor);

    % identify the blocks' class
    [cate_x, cate_y] = block_identify( currentParticle, y_coor, x_coor, fore_x_coor, back_x_coor, fore_y_coor, back_y_coor );

    z_label_cu = zeros(dataNum, dataNum);
    for ii = (1:dataNum)

        % calculate the original coor_xi's value
        forward_sel = (((coor_eta(ii)-fore_y_coor).*(coor_eta(ii)-back_y_coor))<0);     %%% this is for the x coordinate's calculation
        if sum(forward_sel)>0
            % the cut point on the vertical direction
            propo_x_coor = back_x_coor(forward_sel) + ((coor_eta(ii)-back_y_coor(forward_sel)).*(fore_x_coor(forward_sel)-...
                back_x_coor(forward_sel)))./(fore_y_coor(forward_sel)-back_y_coor(forward_sel));
            propo_x_coor = sort(propo_x_coor);
        else
            propo_x_coor = [];
        end

        coor_etaLoc = zeros(1, numel(coor_xi));
        for jj = 1:numel(coor_xi)
            coor_etaLoc(jj) = sum(coor_xi(jj)>(propo_x_coor))+1;
        end

        x_spec = sum(coor_eta(ii)>y_coor)+1;
        origin_tt_cate_x = cate_x(x_spec, coor_etaLoc);

        z_label_cu(ii, :) = origin_tt_cate_x;
    end
    z_label_seq{dd+1} = z_label_cu;
end
