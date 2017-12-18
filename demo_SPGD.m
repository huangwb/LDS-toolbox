clear
close
addPaths; % Add all ducuments to the current path

disp('Set-ups')
feval('param_dyntex_plus');

disp('Preparing dataset')
imdbPath = fullfile(expDir, sprintf('imdb_%s.mat',isplit{1}));
if exist(imdbPath, 'file')
    imdb = load(imdbPath) ;
else
    imdb = setupData(dbDir,imgDir,isplit{1});
    if ~exist(expDir, 'dir')
        mkdir(expDir) ;
    end
    save(imdbPath, '-struct', 'imdb') ;
end

time = zeros(DL_opts.nIter+1,length(DL_opts.nAtoms));
accuracy = zeros(1,length(DL_opts.nAtoms));
costs = zeros(DL_opts.nIter+1,length(DL_opts.nAtoms));
for t=1:length(DL_opts.nAtoms)
    tDL_opts = DL_opts;
    tDL_opts.nAtoms = DL_opts.nAtoms(t);
    
    disp('Learning dictionary atoms -------------------------------------');
%     [dict,this_costs,this_time] = learnDict_Batch(imdb,getBlockBatch,LDS_opts,SC_opts,tDL_opts,STPM);
    [dict,this_costs,this_time] = learnDict_SPGD(imdb,getBlockBatch,LDS_opts,SC_opts,tDL_opts,STPM);
    disp('-----------------------------------------------------------End.')
    time(:,t) = this_time;
    costs(:,t) = this_costs;
    
    % restore K_D to save computations below
    SC_opts.K_D = SC_opts.kernel(dict);
   
    trBatch = find(imdb.images.set==1);
    teBatch = find(imdb.images.set==3);
    % pSTPM.blockPerVid = 100;
    disp('Computing sparse codes for training data-----------------------------')
    trFeatures = zeros(length(dict),length(trBatch));
    for iter=1:length(trBatch)
        thisBlocks = getBlocksPerVid(imdb,trBatch(iter),STPM);
        thisLds = getLdsBatch(thisBlocks.data,LDS_opts);
        trFeatures(:,iter) =  sparse_coding(thisLds,dict,SC_opts);
    end
    disp('-----------------------------------------------------------End.')
    
    disp('Computing sparse codes for testing data-----------------------------')
    teFeatures = zeros(length(dict),length(teBatch));
    for iter=1:length(teBatch)
        thisBlocks = getBlocksPerVid(imdb,teBatch(iter),STPM);
        thisLds = getLdsBatch(thisBlocks.data,LDS_opts);
        teFeatures(:,iter) = sparse_coding(thisLds,dict,SC_opts);
    end
    disp('-----------------------------------------------------------End.')
    
    disp('Training and predicting')
    train_x = trFeatures';
    test_x = teFeatures';
    train_y = imdb.images.label(trBatch)';
    test_y = imdb.images.label(teBatch)';
    [this_accuracy,predict_y] = classify(train_x, train_y, test_x, test_y, Classifier_opts);
    disp('-----------------------------------------------------------End.')
    
    accuracy(t) = this_accuracy;
    SC_opts = rmfield(SC_opts, 'K_D');
end

disp('Saving results.')
% save(fullfile(expDir,[num2str(DL_opts.nAtoms),'_SPGD_Grass-3.mat']),'costs', 'time');