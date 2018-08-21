function currentParticle = Update_particles(TESTFLAG, datas, dataNum, coor_xi, coor_eta, currentParticle, particleN, maxStage, particleBudget, bPara)
% Update the particles


nparticles = Initialize_Particles(particleN, maxStage, particleBudget, dataNum); % initialize the particles
weights = zeros(numel(nparticles), 1);
remainingBudget_copy = ones(numel(nparticles), 1)*particleBudget;


flag_i = 1;
while (max(remainingBudget_copy) > 0)||(flag_i<=numel(currentParticle.ll_ratio_seq))
    
    if flag_i > 1
        resample_likeli = exp(total_weights-max(total_weights));
        % Re-sampling
        sampled_index = randsample(length(resample_likeli),length(resample_likeli)-1, true,transpose(resample_likeli));

        noncurrentParticle_index = find(sampled_index~=1);
        nparticles(noncurrentParticle_index) = nparticles(sampled_index(noncurrentParticle_index)-1);

        nparticles(sampled_index==1) = {truncate_current_particle(currentParticle, flag_i)};

        

        for jj =1:numel(nparticles)
            remainingBudget_copy(jj) = nparticles{jj}.remainingBudget(end);
        end
        
        
    end
    
    for ii = 1:particleN
        if remainingBudget_copy(ii)>0
            kdtree = nparticles{ii}.kdtree;
            nodeNum = nparticles{ii}.nodeNum;
            pointNum = nparticles{ii}.pointNum;
            points = nparticles{ii}.points;
            pIndex = nparticles{ii}.pIndex;

            if length(nparticles{ii}.z_label)<flag_i
                a = 1;
            end
            
            z_label = nparticles{ii}.z_label{flag_i};

            [ kdtree, nodeNum, pointNum, points, pIndex, cut_length, likelihood_ratio, z_label] = propose_cut_in_particle(kdtree, nodeNum, pointNum, points, pIndex, bPara, z_label, datas, coor_xi, coor_eta);
            
            nparticles{ii}.kdtree = kdtree;
            nparticles{ii}.nodeNum = nodeNum;
            nparticles{ii}.pointNum = pointNum;
            nparticles{ii}.points = points;
            nparticles{ii}.pIndex = pIndex;
            nparticles{ii}.z_label{flag_i+1} = z_label;
            
            
            nparticles{ii}.perimeters = nparticles{ii}.perimeters + cut_length;
            nparticles{ii}.ll_ratio_seq = [nparticles{ii}.ll_ratio_seq; likelihood_ratio];
            
            currentCost = -log(1-rand)/(nparticles{ii}.perimeters);
            
            remainingBudget_ii = nparticles{ii}.remainingBudget(end) - currentCost;
            nparticles{ii}.remainingBudget = [nparticles{ii}.remainingBudget; remainingBudget_ii];
            
            % Update the weights for the current nparticle{ii}, the weights' update should simply be the likelihood ratio's update
            weights(ii) = likelihood_ratio;
        else
            weights(ii) =  0;
        end
        
    end
    flag_i = flag_i + 1;
    if flag_i<=numel(currentParticle.ll_ratio_seq)
        fixedWeight = currentParticle.ll_ratio_seq(flag_i);  % fix this and everything will be OK, YEAH!!!
    else
        fixedWeight = 0;
    end
    % Normalize the weights here
    total_weights = [fixedWeight; weights];

end


% Sample one particle for the currentParticle
cP_likili = exp(total_weights-max(total_weights));
cP_index = sum((rand*sum(cP_likili))>[0;cumsum(cP_likili)]);
if cP_index > 1

    if TESTFLAG == 1
        testIndex = currentParticle.testIndex;
        testData = currentParticle.testData;
    end
    currentParticle = nparticles{cP_index-1};
    fprintf('Change the current Particle. \n');
    if TESTFLAG == 1    
        currentParticle.testIndex = testIndex;
        currentParticle.testData = testData;
    end
    
end











end

