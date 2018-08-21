function [currentParticle, datas] = Initialize_currentParticle(TESTFLAG, dataNum, particleBudget, maxStage, bPara, datas, coor_xi, coor_eta)
% Initialize one particle
% dataNum: the number of nodes in the network
% particleBudget: the budget parameter for the BSP-Tree Process
% maxStage: maximum number of stages for PG in inferring the BSP-Tree
% Process relational model
% bPara: the parameters of Beta distribution in the likelihood function

kdtree = zeros(maxStage, 4); % define the kd-tree structure for each particle
kdtree(1, :) = [1 0 0 4]; % kd-tree column meaning: 1, node ID; 2, 3, child IDs; 4, perimeter of the node

nodeNum = 1;
pointNum = 4;
points = zeros((maxStage+4), 2);
points(1:4, :) = [0, 0;1, 0; 1, 1; 0, 1];
pIndex = cell(maxStage, 1);
pIndex{1} = 1:4;


if TESTFLAG==1
    testMat = zeros(dataNum, dataNum);
    for jj = 1:dataNum
        testMat(jj, randsample(dataNum, floor(dataNum/10))) = 1;
    end
    test_data = datas(testMat==1);
    datas(testMat==1) = -1;
end


remainingBudget = particleBudget;
currentPerimeter = 4;
ll_ratio_seq = 0;
remainingBudget_seq = particleBudget;
z_label = ones(size(datas));

z_label_seq = cell(500, 1);

flag_i = 1;
while remainingBudget > 0
    
    % select the polygon to be cut

    [ kdtree, nodeNum, pointNum, points, pIndex, cut_length, likelihood_ratio, z_label] = propose_cut_in_particle(kdtree, nodeNum, pointNum, points, pIndex, bPara, z_label, datas, coor_xi, coor_eta);    

    % need to use the data as the likelihood to update weight
    ll_ratio_seq = [ll_ratio_seq; likelihood_ratio];
    
    % generate_cost from the Exponential distribution
    currentPerimeter = currentPerimeter + 2*cut_length;
    currentCost = -log(1-rand)/(currentPerimeter);
    
    remainingBudget = remainingBudget - currentCost;
    remainingBudget_seq = [remainingBudget_seq; remainingBudget];
    z_label_seq{flag_i} = z_label;
    flag_i = flag_i + 1;
    
end

currentParticle.kdtree = kdtree;
currentParticle.nodeNum = nodeNum;
currentParticle.points = points;
currentParticle.pIndex = pIndex;
currentParticle.pointNum = pointNum;

currentParticle.ll_ratio_seq = ll_ratio_seq;
currentParticle.perimeters = currentPerimeter;
currentParticle.remainingBudget = remainingBudget_seq;

currentParticle.z_label = z_label_seq;

currentParticle.testIndex = find(testMat==1);
currentParticle.testData = test_data;



end

