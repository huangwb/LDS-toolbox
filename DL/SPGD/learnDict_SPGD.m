function [dict,costs,time] = learnDict_SPGD(imdb,getBlockBatch,LDS_opts,SC_opts,DL_opts)
%SPGD performs dictionary learning on LDSs with the Stochastic-Projected-Gradient-Decent method
%
% INUTS
% imdb            - {} dataset
% getBatch   - The function handle to get training batch from imdb
% DL_opts      - . Dictionary learning set-up parameters, see ./params/param_default.m
% SC_opts      - . Sparse coding set-up parameters, see ./params/param_default.m
%
% OUTPUTS
% dict        - {} Learned dictionary
% implemented by Wenbing Huanng, 2016-12-14

rng('default');
rng(0) ;

% set training, testing and validation sets
trBatch = find(imdb.images.set==1);
teBatch = find(imdb.images.set==3);
% tp = randperm(size(teBatch,2));
% teBatch = teBatch(tp(1:DL_opts.teBatchSize));

disp('Setting validation batches.');
teBlocks = getBlockBatch(imdb,teBatch);
teLdsBatch = getLdsBatch(teBlocks.data,LDS_opts);

% begin dictionary learning
rng('default');
rng(0) ;
% initialize the dictionary
switch DL_opts.initial_method
    case 'RAND'
        disp('Initializing dictionary atoms with random video blocks.');
        initBatch = trBatch(randperm(length(trBatch)));
        initBlocks = getBlockBatch(imdb,initBatch);
        initLds = getLdsBatch(initBlocks.data,LDS_opts);
%         dict = initLds(1:DL_opts.nAtoms);
        rand_idx = randperm(length(initLds));
        dict = initLds(rand_idx(1:DL_opts.nAtoms));
    case 'KMEAN'
        fprintf('Initializing the dictionary using kmeans algorithm.\n');
%         centers = kmeans_LDS_SymSkew(data,RandDict,DL_opts.nAtoms,LDS_opts,5);
%         InitDict = cell(1,2*DL_opts.nAtoms);
%         for j=1:DL_opts.nAtoms
%             InitDict{j} = centers{j};
%         end
end

switch LDS_opts.decomposition
    case 'CANONICAL'
        update_dict = @update_Can;
    case 'SYMSKEW'
        update_dict = @update_SymSkew;
end

SPGD_opts.SN_a = LDS_opts.SN_a;
SPGD_opts.SN_b = LDS_opts.SN_b;
SPGD_opts.learningRate_A = DL_opts.learningRate_A;
SPGD_opts.learningRate_C = DL_opts.learningRate_C;

% compute the initial cost for testing batch
SC_opts.L = min(int32(DL_opts.nAtoms/2),20);
costs = zeros(DL_opts.nIter+1,1);
time = zeros(DL_opts.nIter+1,1);
t1=clock;
teAlpha = sparse_coding(teLdsBatch,dict,SC_opts);
% teAlpha = local_sparse_coding(teLdsBatch,dict,SC_opts);
costs(1) = computeDictCost(teLdsBatch,dict,teAlpha,SC_opts);
t2 = clock;
t = etime(t2,t1);
fprintf('Initial objective cost -->%.3f; time cost -->%.3fs\n',costs(1),t);

for iter = 1:DL_opts.nIter
    t1=clock;
    nVideo = length(trBatch);
    vBatch = trBatch(randperm(nVideo));
    
    for j=1:ceil(nVideo/DL_opts.vBatchSize)
        tj1 = clock;
        vBatch_start = (j-1)*DL_opts.vBatchSize+1;
        vBatch_end = min(vBatch_start + DL_opts.vBatchSize - 1,nVideo);
        blocks = getBlockBatch(imdb,vBatch(vBatch_start:vBatch_end));
        ldsBatch = getLdsBatch(blocks.data,LDS_opts);
        
        % update codes and dictonary atoms
        trAlpha = sparse_coding(ldsBatch,dict,SC_opts);
%         trAlpha = local_sparse_coding(ldsBatch,dict,SC_opts);
        dict = update_dict(ldsBatch,dict,trAlpha,SPGD_opts);
        tj2 = clock;
        tj = etime(tj2,tj1);
        fprintf('Iter#%d: %d/%d, time cost -->%.3fs\n',iter, j, ceil(nVideo/DL_opts.vBatchSize), tj);
        
    end
    
    % compute the testing objective after each update
    teAlpha = sparse_coding(teLdsBatch,dict,SC_opts);
%     teAlpha = local_sparse_coding(teLdsBatch,dict,SC_opts);
    costs(iter+1) = computeDictCost(teLdsBatch,dict,teAlpha,SC_opts);
    
    t2 = clock;
    t = etime(t2,t1);
    time(iter+1) = t;
    fprintf('Iter#%d: objective cost -->%.3f; time cost -->%.3fs\n',iter,costs(iter+1),t);
    
    if mod(iter, DL_opts.learningRate_stepsize) == 0
        SPGD_opts.learningRate_A = DL_opts.learningRate_gamma * SPGD_opts.learningRate_A;
        SPGD_opts.learningRate_C = DL_opts.learningRate_gamma * SPGD_opts.learningRate_C;
    end
    
end

fprintf('-------\n');

% save results
% ...
end