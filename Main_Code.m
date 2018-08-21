% Main Code for Binary Space Partitioning-Tree Process Relational Model


% clear;
% clc;

% load data: datas, dataNum

load('flickr_subset.mat');

% indicate if the data is split into training data and test data
TESTFLAG = 1;

total_iter = 10;

% parameters for the Global MCMC setting
maxIter = 20;

% start to initialize particleN particles
% parameters for the Particles
particleN = 10;
maxStage = 500;
particleBudget = 0.8;
bPara = [1, 1];
internalRound = 10;

% initialzie the currentParticle and model parameters
coor_xi = rand(dataNum, 1);
coor_eta = rand(dataNum, 1);

[currentParticle, datas] = Initialize_currentParticle(TESTFLAG, dataNum, particleBudget, maxStage, bPara, datas, coor_xi, coor_eta);


trainLLSEQ = zeros(maxIter, 1);
train_auc_SEQ = zeros(maxIter, 1);

for i1 = 1:maxIter
    
    % Update nodes' coordinates
    for j1 =1:internalRound
        [coor_xi, coor_eta, currentParticle] = Update_coordinates(TESTFLAG, currentParticle, coor_xi, coor_eta, datas, dataNum, bPara);    
    end
    trainLLSEQ(i1) = currentParticle.trainll;
    train_auc_SEQ(i1) = currentParticle.aucValue;
    train_auc_SEQ(i1)
    % Update nparticles
    currentParticle = Update_particles(TESTFLAG, datas, dataNum, coor_xi, coor_eta, currentParticle, particleN, maxStage, particleBudget, bPara);    
    
end


figure(1);
h4 = plot(trainLLSEQ(1:(i1)));
figure(2);
h5 = plot(train_auc_SEQ(1:(i1-1)));

