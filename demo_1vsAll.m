clear 
close 
addPaths; % Add all ducuments to the current path

disp('Set-ups')
feval('param_dyntex11');

accuracy_dict = zeros(numel(LDS_opts.dim_hidden),numel(isplit));
time_dict = zeros(numel(LDS_opts.dim_hidden),numel(isplit));
%load temporal results
% load(fullfile(expDir,[num2str(DL_opts.nAtoms),'.mat']));
for n=1:numel(LDS_opts.dim_hidden)
    
    thisLDS_opts = LDS_opts;
    thisLDS_opts.dim_hidden = LDS_opts.dim_hidden(n);
    imdbPath = fullfile(expDir, sprintf('imdb_%s.mat',isplit{1}));
    t1=clock;
    disp('Preparing dataset')
    if exist(imdbPath, 'file')
        imdb = load(imdbPath) ;
    else
        imdb = setupData(dbDir,imgDir,isplit);
        if ~exist(expDir, 'dir')
            mkdir(expDir) ;
        end
        save(imdbPath, '-struct', 'imdb') ;
    end
    
    
    disp('Dictionary learning--------------------------------------------');
%     dict = learnDict_Batch(imdb,getBlockBatch,LDS_opts,SC_opts,DL_opts,STPM);
    dict = learnDict_SPGD(imdb,getBlockBatch,LDS_opts,SC_opts,DL_opts,STPM);
    disp('-----------------------------------------------------------End.')
    
    % restore K_D to save computations below
    SC_opts.K_D = SC_opts.kernel(dict);
    
    
    batch = find(imdb.images.set==1);
    pSTPM = STPM;
    % pSTPM.blockPerVid = 100;
    disp('Pooling features for all data-----------------------------')
    features = zeros(length(dict)*STPM.nBins,length(batch));
    for iter=1:length(batch)
        tic;
        thisBlocks = getBlocksPerVid(imdb,batch(iter),STPM);
        %     thisBlocks = getBlockBatch(imdb,trBatch(iter),pSTPM);
        features(:,iter) = SCPooling(thisBlocks, dict, pSTPM.pyramid, SC_opts, thisLDS_opts);
        t = toc;
        fprintf('Processing training data %d takes %f seconds\n',iter, t);
    end
    disp('-----------------------------------------------------------End.')
    
    nClasses = length(imdb.classes.name);
    nData = length(batch);
    outputs = zeros(nClasses,nData);
    for f=1:length(batch)
        teBatch = f;
        teFeatures = features(:,teBatch);
        trBatch = setdiff(1:length(batch),f);
        trFeatures = features(:, trBatch);
        
        disp('Training and predicting')
        [accuracy,predict_y] = classify(trFeatures', imdb.images.label(trBatch)', teFeatures', imdb.images.label(teBatch)', Classifier_opts);
        disp('-----------------------------------------------------------End.')
    
        outputs(predict_y,f) = 1;
        accuracy_dict(n,f) = accuracy;
        t2=clock;
        time_dict(n,f) = etime(t2,t1);
        
%         fprintf('Fold %d: Accuracy with a learned dictionary of size %d with hidden order %d: %.1f%%; Time cost: %f\n',f, DL_opts.nAtoms,thisLDS_opts.dim_hidden,accuracy_dict(n,f),time_dict(n,f));
        

    end
    mean(accuracy_dict)
    % save results
    disp('Saving results.')
    save(fullfile(expDir,[num2str(DL_opts.nAtoms),'_STPM_can.mat']),'accuracy_dict','outputs', 'features', 'dict', 'LDS_opts','SC_opts','DL_opts');
end


%exit;