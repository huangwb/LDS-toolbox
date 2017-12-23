    function [centers, cluster_idx] = Kmeans(data,opts)
%COMPUTECODEBOOK performs k-means
%
% INUTS
% data            - {} LDS data
%
% OUTPUTS
% centers        - {} Resulting clustering centers
% implemented by Wenbing Huanng, 2016-12-19

rng('default');
rng(0);

MinCostVariation = 1e-4;

nData = length(data);
computeKernel = opts.kernel;
computeMean = opts.computeMean;

mean_opts = opts;
mean_opts.nIter = 5;

cluster_idx = randi(opts.nAtoms,[1,nData]);
centers = cell(1,opts.nAtoms);
for iter = 1:opts.nIter
    fprintf('.');
    t1=clock;
    for i = 1:opts.nAtoms
        idx = find(cluster_idx == i);
        if (isempty(idx))
            %zombie centers
            randVal = randperm(nData);
            centers{i} = data{randVal(1)};
        elseif iter ==1
            centers{i} = computeMean(data(idx),mean_opts);
		elseif iter > 1
			centers{i} = computeMean(data(idx),mean_opts, centers{i});
        end
    end
    
    %assign points and compute the cost
    [currCost,cluster_idx] = kmeans_cost(data,centers,computeKernel);
    
    t2=clock;
    time_cost = etime(t2,t1);
    
    if (iter == 1)
        fprintf('kmeans: initial cost is %6.3f, time cost %6.3f.\n',currCost, time_cost);
    else
        cost_diff = norm(preCost - currCost) ;
        if (cost_diff < MinCostVariation)
            break ;
        else
            fprintf('kmeans: Iter#%d, cost is %6.3f, time cost is %6.3f.\n',iter,currCost, time_cost);
        end
    end
    preCost = currCost;
    
    % decrease the learning rate
%     if mod(iter, opts.learningRate_stepsize) == 0
%         mean_opts.learningRate_A = mean_opts.learningRate_gamma * mean_opts.learningRate_A;
%         mean_opts.learningRate_C = mean_opts.learningRate_gamma * mean_opts.learningRate_C;
%     end
    
end
fprintf('\n');
end

function [outCost,cluster_idx] = kmeans_cost(data,centers,computeKernel)

nData = length(data);
l_sim = computeKernel(centers,data);
[maxSim,cluster_idx] = max(l_sim);
outCost = 1/nData*sum(maxSim);

end

