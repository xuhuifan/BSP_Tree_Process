function [ kdtree, nodeNum, pointNum, points, pIndex, cut_length, likelihood_ratio, z_label] = propose_cut_in_particle( kdtree, nodeNum, pointNum, points, pIndex, bPara, z_label, datas, coor_xi, coor_eta)
% Propose a cut upon the current partitioning structure ()

cutCandidate = find((kdtree(:, 1)>0)&(kdtree(:, 2)==0));
leafNum = numel(cutCandidate);

polygonPerimeter = kdtree(cutCandidate, 4);
cutIndex = sum(rand*sum(polygonPerimeter)>[0;cumsum(polygonPerimeter)]);

polygonIndex = kdtree(cutCandidate(cutIndex), 1);

sequentialPointsIndex = pIndex{polygonIndex};
sequentialPoints = points(sequentialPointsIndex, :);


%%% propose to cut the Polygon with the current method
%(the angle are assumed to be uniformly
% sampled, which is not consistent with the generative model. However, since this is only an
% initialization, this should be OK. ANYWAY, EXTRA CARE SHOULD BE TAKEN IN THIS CASE. )

    [intersects_two, intersects_index] = sequentialPropose( sequentialPoints );

if intersects_index(2)<numel(sequentialPointsIndex)
    pIndex_new1 = [sequentialPointsIndex(1:intersects_index(1)), pointNum+1, pointNum+2, sequentialPointsIndex((intersects_index(2)+1):end)];
    pIndex_new2 = [pointNum+1, sequentialPointsIndex((intersects_index(1)+1):intersects_index(2)), pointNum+2];
else
    pIndex_new1 = [sequentialPointsIndex(1:intersects_index(1)), pointNum+1, pointNum+2];
    pIndex_new2 = [pointNum+1, sequentialPointsIndex((intersects_index(1)+1):intersects_index(2)), pointNum+2];
end

% update kdtree structure
points(pointNum+1, :) = intersects_two(1, :);
points(pointNum+2, :) = intersects_two(2, :);
pIndex{nodeNum+1} = pIndex_new1;
pIndex{nodeNum+2} = pIndex_new2;

kdtree(polygonIndex, 2) = nodeNum+1;
kdtree(polygonIndex, 3) = nodeNum+2;

kdtree(nodeNum+1, 1) = nodeNum+1;
kdtree(nodeNum+1, 4) = perimeter_cal(points(pIndex_new1, :));
kdtree(nodeNum+2, 1) = nodeNum+2;
kdtree(nodeNum+2, 4) = perimeter_cal(points(pIndex_new2, :));

nodeNum = nodeNum + 2;
pointNum = pointNum + 2;

cut_length = sqrt(sum((intersects_two(1, :) - intersects_two(2, :)).^2));

% start to calculate the ratio between the new likelihood and the old
% likelihood

polygon_1_points = points(pIndex_new1, :); % the points along the fist new polygon

k_index = find(z_label==cutIndex); % check out the nodes' with the selected polygon

xis = repmat(coor_xi', size(datas, 1), 1);
etas = repmat(coor_eta, 1, size(datas, 1));
xi_index = xis(k_index);
eta_index = etas(k_index);


prop_z = z_label;
index_points = [xi_index,eta_index];
k1_set = zeros(1, numel(k_index));
for ii = 1:numel(k_index) % to check the remaining data's belonging one by one
    kk_origin = bsxfun(@minus, polygon_1_points, index_points(ii, :));
    
    if all((kk_origin([2:end 1], 1).*kk_origin(:, 2)-kk_origin(:, 1).*kk_origin([2:end 1], 2))<0)  % maybe this is right, let us see
        k1_set(ii) = ii;
    end
end

k1_set(k1_set==0) = [];
k1_point = k_index(k1_set);
prop_z(k1_point) = leafNum+1;

k_index(k1_set) = [];
k2_point = k_index;
prop_z(k2_point) = leafNum+2;

%%%
%%%

point1 = datas(k1_point);
point2 = datas(k2_point);

tau_kl1 = [sum(sum(point1==1)) sum(sum(point1==0))];
tau_kl2 = [sum(sum(point2==1)) sum(sum(point2==0))];
tau_kl0 = tau_kl1 + tau_kl2;

likelihood_ratio = (-(sum(gammaln(tau_kl0+bPara))-gammaln(sum(bPara)+sum(tau_kl0))+...
    sum(gammaln(bPara))-gammaln(sum(bPara))+...
    gammaln(sum(bPara)+sum(tau_kl1))-sum(gammaln(tau_kl1+bPara))+...
    gammaln(sum(bPara)+sum(tau_kl2))-sum(gammaln(tau_kl2+bPara))));  %%%%%%% to be simplified later

z_label = prop_z;
z_label(z_label>cutIndex) = z_label(z_label>cutIndex) - 1;

end

