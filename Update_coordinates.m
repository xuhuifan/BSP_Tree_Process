function [coor_xi, coor_eta, currentParticle] = Update_coordinates(TESTFLAG, currentParticle, coor_xi, coor_eta, datas, dataNum, bPara)
% Update the coordinates
% TESTFLAG: indicate if we need to calculate the AUC and test log-likelihood
% currentParticle: the current structure of particle
% coor_xi, coor_eta: row and column coordinates for all the nodes
% datas, dataNum, bPara: defined as previously


% use Metropolis-Hasting Algorithm to update the nodes' coordinates

points = currentParticle.points;
kdtree = currentParticle.kdtree;
nodeNum = currentParticle.nodeNum;


%
%
% model paramters

z_label = currentParticle.z_label{(nodeNum-1)/2};

%
%
%
numClass = (nodeNum+1)/2;

tau1_kl = zeros(1, numClass);
tau0_kl = zeros(1, numClass);
for ii = 1:numClass
    ii_data = datas(z_label == ii);
    tau1_kl(ii) = sum(sum(ii_data==1));
    tau0_kl(ii) = sum(sum(ii_data==0));
end


%% First, we need to check out the available points
if nodeNum > 1
    
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
    
    
    for ii = randperm(dataNum)
        %& First sampling \coor_xi_ii
        propo_coor_eta = rand;

        
        origin_ii1 = hist((z_label(ii, datas(ii, :)==1)), 1:numClass);
        origin_ii0 = hist((z_label(ii, datas(ii, :)==0)), 1:numClass);
        
        
        tau1_kl_ori = tau1_kl - origin_ii1;
        tau0_kl_ori = tau0_kl - origin_ii0;
        
        
        % calculate the proposal_coor_xi's value
        forward_sel = (((propo_coor_eta-fore_y_coor).*(propo_coor_eta-back_y_coor))<0);     %%% this is for the x coordinate's calculation
        if sum(forward_sel)>0
            % the cut point on the vertical direction
            propo_x_coor = back_x_coor(forward_sel) + ((propo_coor_eta-back_y_coor(forward_sel)).*(fore_x_coor(forward_sel)-...
                back_x_coor(forward_sel)))./(fore_y_coor(forward_sel)-back_y_coor(forward_sel));
            propo_x_coor = sort(propo_x_coor);
        else
            propo_x_coor = [];
        end
        
        coor_etaLoc = zeros(1, numel(coor_xi));
        for jj = 1:numel(coor_xi)
            coor_etaLoc(jj) = sum(coor_xi(jj)>(propo_x_coor))+1;
        end
        
        x_spec = sum(propo_coor_eta>y_coor)+1;
        tt_cate_x = cate_x(x_spec, coor_etaLoc);
        
        % start to define the original categories
        
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
        %
        
        propo_ii1 = hist((tt_cate_x.*(datas(ii, :)==1)), 0:numClass);
        propo_ii1(1) = [];
        propo_ii0 = hist((tt_cate_x.*(datas(ii, :)==0)), 0:numClass);
        propo_ii0(1) = [];
        tau1_propo = tau1_kl_ori + propo_ii1;
        tau0_propo = tau0_kl_ori + propo_ii0;
        
        if (sum(tau1_propo<0)>0)||(sum(tau0_propo<0)>0)
            fprintf('xi_sampling wrong, propro\n');
        end
        
        
        
        
        if log(rand)<(sum(gammaln(bPara(1)+tau1_propo)+gammaln(bPara(2)+...
                tau0_propo)-gammaln(bPara(1)+tau1_propo+bPara(2)+tau0_propo)-...
                gammaln(bPara(1)+tau1_kl)-gammaln(bPara(2)+tau0_kl)+...
                gammaln(bPara(1)+tau1_kl+bPara(2)+tau0_kl)))
            tau1_kl = tau1_propo;
            tau0_kl = tau0_propo;
            coor_eta(ii) = propo_coor_eta;
            z_label(ii, :) = tt_cate_x;
            
        else
            z_label(ii, :) = origin_tt_cate_x;
            
            propo_ii1 = hist((origin_tt_cate_x.*(datas(ii, :)==1)), 0:numClass);
            propo_ii1(1) = [];
            propo_ii0 = hist((origin_tt_cate_x.*(datas(ii, :)==0)), 0:numClass);
            propo_ii0(1) = [];
            tau1_kl = tau1_kl_ori + propo_ii1;
            tau0_kl = tau0_kl_ori + propo_ii0;
        end
        
        
        
        %% Then sampling \coor_eta_ii
        propo_coor_xi = rand;
        
        % delete the current statistics from Bk1, Bk0
        
        origin_ii1 = hist((z_label(datas(:, ii)==1, ii)), 1:numClass);
        origin_ii0 = hist((z_label(datas(:, ii)==0, ii)), 1:numClass);
        
        tau1_kl_ori = tau1_kl - origin_ii1;
        tau0_kl_ori = tau0_kl - origin_ii0;
        
        % calculate the proposal_coor_eta's value
        backward_sel = ((propo_coor_xi-fore_x_coor).*((propo_coor_xi-back_x_coor))<0);        %%% this is for the y coordinate's calculation
        
        if sum(backward_sel)>0
            % the cut point on the horizonal direction
            propo_y_coor = back_y_coor(backward_sel) + ((propo_coor_xi-back_x_coor(backward_sel)).*(fore_y_coor(backward_sel)...
                -back_y_coor(backward_sel)))./(fore_x_coor(backward_sel)-back_x_coor(backward_sel));
            propo_y_coor = sort(propo_y_coor);
        else
            propo_y_coor = [];
        end
        
        coor_xiLoc = zeros(1, numel(coor_eta));
        for jj = 1:numel(coor_eta)
            coor_xiLoc(jj) = sum(coor_eta(jj)>(propo_y_coor))+1;
        end
        
        x_spec = sum(propo_coor_xi>x_coor)+1;
        tt_cate_y = cate_y(x_spec, coor_xiLoc);
        
        
        % start to define the original column
        backward_sel = ((coor_xi(ii)-fore_x_coor).*((coor_xi(ii)-back_x_coor))<0);        %%% this is for the y coordinate's calculation
        
        if sum(backward_sel)>0
            % the cut point on the horizonal direction
            propo_y_coor = back_y_coor(backward_sel) + ((coor_xi(ii)-back_x_coor(backward_sel)).*(fore_y_coor(backward_sel)...
                -back_y_coor(backward_sel)))./(fore_x_coor(backward_sel)-back_x_coor(backward_sel));
            propo_y_coor = sort(propo_y_coor);
        else
            propo_y_coor = [];
        end
        
        coor_xiLoc = zeros(1, numel(coor_eta));
        for jj = 1:numel(coor_eta)
            coor_xiLoc(jj) = sum(coor_eta(jj)>(propo_y_coor))+1;
        end
        
        x_spec = sum(coor_xi(ii)>x_coor)+1;
        origin_tt_cate_y = cate_y(x_spec, coor_xiLoc);
        %
        
        propo_ii1 = hist((tt_cate_y'.*(datas(:, ii)==1)), 0:numClass);
        propo_ii0 = hist((tt_cate_y'.*(datas(:, ii)==0)), 0:numClass);
        propo_ii1(1) = [];
        propo_ii0(1) = [];
        tau1_propo = tau1_kl_ori + propo_ii1;
        tau0_propo = tau0_kl_ori + propo_ii0;
        
        if (sum(tau1_propo<0)>0)||(sum(tau0_propo<0)>0)
            fprintf('eta_sampling wrong, propo\n');
        end
        
    
        
        if log(rand)<(sum(gammaln(bPara(1)+tau1_propo)+gammaln(bPara(2)+...
                tau0_propo)-gammaln(bPara(1)+tau1_propo+bPara(2)+tau0_propo)-...
                gammaln(bPara(1)+tau1_kl)-gammaln(bPara(2)+tau0_kl)+...
                gammaln(bPara(1)+tau1_kl+bPara(2)+tau0_kl)))
            
            tau1_kl = tau1_propo;
            tau0_kl = tau0_propo;
            
            coor_xi(ii) = propo_coor_xi;
            z_label(:, ii) = tt_cate_y';
            
            
        else
            z_label(:, ii) = origin_tt_cate_y';
            
            propo_ii1 = hist((origin_tt_cate_y'.*(datas(:, ii)==1)), 0:numClass);
            propo_ii0 = hist((origin_tt_cate_y'.*(datas(:, ii)==0)), 0:numClass);
            propo_ii1(1) = [];
            propo_ii0(1) = [];
            tau1_kl = tau1_kl_ori + propo_ii1;
            tau0_kl = tau0_kl_ori + propo_ii0;           
        end
        
        
    end
    
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %% Start to calculate the AUC value, which is critical
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    
    if TESTFLAG==1
        tau_kl = tau1_kl + tau0_kl;
        linkProb = (tau1_kl+bPara(1))./(tau_kl+sum(bPara));
        
        testIndex = currentParticle.testIndex;
        testData = currentParticle.testData;
        
        predictProb = zeros(1, numel(testIndex));
        
        
        for aucii = 1:numel(testIndex)
            iiRow = mod(testIndex(aucii), dataNum);
            if (iiRow==0)
                iiRow = dataNum;
            end
            iiColumn = ceil(testIndex(aucii)/dataNum);
            
            ii_coor_eta = coor_eta(iiRow);
            ii_coor_xi = coor_xi(iiColumn);
            
            
            %%% From the tt_cate_y perspective
            
            % calculate the proposal_coor_eta's value
            backward_sel = ((ii_coor_xi-fore_x_coor).*((ii_coor_xi-back_x_coor))<0);        %%% this is for the y coordinate's calculation
            
            if sum(backward_sel)>0
                % the cut point on the horizonal direction
                propo_y_coor = back_y_coor(backward_sel) + ((ii_coor_xi-back_x_coor(backward_sel)).*(fore_y_coor(backward_sel)...
                    -back_y_coor(backward_sel)))./(fore_x_coor(backward_sel)-back_x_coor(backward_sel));
                propo_y_coor = sort(propo_y_coor);
            else
                propo_y_coor = [];
            end
            
            coor_xiLoc = sum(ii_coor_eta>(propo_y_coor))+1;
            x_spec = sum(ii_coor_xi>x_coor)+1;
            tt_cate_y = cate_y(x_spec, coor_xiLoc);
            
            %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%
            
            predictProb(aucii) = linkProb(tt_cate_y);
        end        
        
        n1 = sum(testData(:));
        no = numel(testData)-n1;
        
        [~, rank_indcs] = sort(predictProb);
        R_sorted = testData(rank_indcs);
        So = sum(find(R_sorted > 0));
        if ((So - (n1*(n1+1))/2)<0)||(n1<0)||(no<0)
            a = 1;
            fprintf('error. \n');
        end
        aucValue = (So - (n1*(n1+1))/2)/(n1*no);
    else
        aucValue = 0.5;
    end
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %% End the auc value calculation here
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%    %%%%%%%%
    
else
    
    z_label = ones(dataNum, dataNum);
    aucValue = 0.5;
end


tau_kl = tau1_kl + tau0_kl;
probs = (tau1_kl+bPara(1))./(tau_kl+sum(bPara));

trainll = 0;

for ii = 1:numClass
    ii_data = datas(z_label == ii);
    trainll = trainll + sum(ii_data==1)*log(probs(ii))+sum(ii_data==0)*log(1-probs(ii));
end



%
currentParticle.tau1_kl = tau1_kl;
currentParticle.tau0_kl = tau0_kl;
currentParticle.trainll = trainll;
currentParticle.aucValue = aucValue;

currentParticle.z_label = z_label_determine(currentParticle, coor_eta, coor_xi, points, kdtree, nodeNum, dataNum);

end

