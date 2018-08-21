function nparticles = Initialize_Particles(particleN, maxStage, particleBudget, dataNum)
% Initialize the particleN particles

nparticles = cell(particleN, 1);

points = zeros((maxStage+4), 2);
points(1:4, :) = [0, 0;1, 0; 1, 1; 0, 1];
pIndex = cell(maxStage, 1);
pIndex{1} = 1:4;
pointNum = 4;

for ii = 1:particleN
    nparticles{ii}.kdtree = zeros(maxStage, 4); % define the kd-tree structure for each particle
    nparticles{ii}.kdtree(1, :) = [1 0 0 4]; % kd-tree column meaning: 1, node ID; 2, 3, child IDs; 4, perimeter of the node
    
    nparticles{ii}.nodeNum = 1;
    nparticles{ii}.pointNum = pointNum;
    
    nparticles{ii}.points = points;
    nparticles{ii}.pIndex = pIndex;
    
    nparticles{ii}.budget = particleBudget;
    nparticles{ii}.cost = 0;
    
    nparticles{ii}.ll_ratio_seq = 1;
    nparticles{ii}.perimeters = 4;
    
    nparticles{ii}.z_label = {ones(dataNum, dataNum)};
    
    nparticles{ii}.remainingBudget = particleBudget;
    
end



end

