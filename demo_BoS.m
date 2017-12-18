clear
close
addPaths; % Add all ducuments to the current path

rng('default');
rng(0);

disp('Set-ups')
feval('param_dyntex');

for n=1:numel(LDS_opts.dim_hidden)
    t1=clock;
    disp('Preparing dataset')
    imdbPath = fullfile(expDir, sprintf('imdb_%s.mat',num2str(isplit{1})));
    if exist(imdbPath, 'file')
        imdb = load(imdbPath) ;
    else
        imdb = setupData(dbDir,imgDir,isplit{f});
        if ~exist(expDir, 'dir')
            mkdir(expDir) ;
        end
        save(imdbPath, '-struct', 'imdb') ;
    end
    trBatch = find(imdb.images.set==1);
    teBatch = find(imdb.images.set==3);
    
    disp('Clustering via K-means-----------------------------')
    trBlocks = getBlockBatch(imdb,trBatch,BoS_opts);
    trainLds = getLdsBatch(trBlocks.data,LDS_opts);
    centers = Kmeans(trainLds, BoS_opts);
    disp('-----------------------------------------------End.')
    
%     batch = find(imdb.images.set==1);
%     % pSTPM.blockPerVid = 100;
%     disp('Pooling features for all data-----------------------------')
%     features = zeros(length(batch), length(centers));
%     for iter=1:length(batch)
%         tic;
%         thisBlocks = getBlocksPerVid(imdb,batch(iter),BoS_opts);
%         thisLds = getLdsBatch(thisBlocks.data,LDS_opts);
%         features(iter,:) = BoSCoding(thisLds, centers, BoS_opts);
%         t = toc;
%         fprintf('Processing training data %d takes %f seconds\n',iter, t);
%     end
%     disp('-----------------------------------------------------------End.')
%     
%     accuracy_dict = zeros(numel(LDS_opts.dim_hidden),numel(isplit));
%     time_dict = zeros(numel(LDS_opts.dim_hidden),numel(isplit));
%     nClasses = length(imdb.classes.name);
%     nData = length(batch);
%     outputs = zeros(nClasses,nData);
%     for f=1:length(batch)
%         teBatch = f;
%         teFeatures = features(teBatch,:);
%         trBatch = setdiff(1:length(batch),f);
%         trFeatures = features(trBatch,:);
%         
%         disp('Training and predicting')
%         [accuracy,predict_y] = classify(trFeatures, imdb.images.label(trBatch)', teFeatures, imdb.images.label(teBatch)', Classifier_opts);
%         disp('-----------------------------------------------------------End.')
%         
%         outputs(predict_y,f) = 1;
%         accuracy_dict(n,f) = accuracy;
%         t2=clock;
%         time_dict(n,f) = etime(t2,t1);
%         
%         %         fprintf('Fold %d: Accuracy with a learned dictionary of size %d with hidden order %d: %.1f%%; Time cost: %f\n',f, DL_opts.nAtoms,thisLDS_opts.dim_hidden,accuracy_dict(n,f),time_dict(n,f));
%         
%         
%     end
%     mean(accuracy_dict)

disp('Computing codes for training data-----------------------------')
trFeatures = zeros(length(trBatch),BoS_opts.nAtoms);
for iter=1:length(trBatch)
    tic;
    thisBlocks = getBlocksPerVid(imdb,trBatch(iter),BoS_opts);
    thisLds = getLdsBatch(thisBlocks.data,LDS_opts);
    trFeatures(iter,:) = BoSCoding(thisLds, centers, BoS_opts);
    t = toc;
    fprintf('Processing training data %d takes %f seconds\n',iter, t);
end
disp('-----------------------------------------------------------End.')

disp('Computing codes for testing data-----------------------------')
teFeatures = zeros(length(teBatch),BoS_opts.nAtoms);
for iter=1:length(teBatch)
    tic;
    thisBlocks = getBlocksPerVid(imdb,teBatch(iter),BoS_opts);
    thisLds = getLdsBatch(thisBlocks.data,LDS_opts);
    teFeatures(iter,:) = BoSCoding(thisLds, centers, BoS_opts);
    t = toc;
    fprintf('Processing testing data %d takes %f seconds\n',iter, t);
end
disp('-----------------------------------------------------------End.')

disp('Training and predicting')
train_x = trFeatures;
test_x = teFeatures;
train_y = imdb.images.label(trBatch)';
test_y = imdb.images.label(teBatch)';
[accuracy,predict_y] = classify(train_x, train_y, test_x, test_y, Classifier_opts);
accuracy
disp('-----------------------------------------------------------End.')

end

disp('Saving results.')
save(fullfile(expDir,[num2str(DL_opts.nAtoms),'_BoS_Grass.mat']),'accuracy', 'centers', 'LDS_opts','SC_opts','BoS_opts');
