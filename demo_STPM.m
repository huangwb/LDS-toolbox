clear
close
addPaths; % Add all ducuments to the current path

disp('Set-ups')
feval('param_dyntex');

t1=clock;

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

disp('Learning dictionary atoms -------------------------------------');
%dict = learnDict_Batch(imdb,getBlockBatch,LDS_opts,SC_opts,DL_opts,STPM);
dict = learnDict_SPGD(imdb,getBlockBatch,LDS_opts,SC_opts,DL_opts,STPM);
disp('-----------------------------------------------------------End.')

% restore K_D to save computations below
SC_opts.K_D = SC_opts.kernel(dict);

trBatch = find(imdb.images.set==1);
teBatch = find(imdb.images.set==3);

% pSTPM.blockPerVid = 100;
disp('Pooling features for training data-----------------------------')
trFeatures = zeros(length(dict)*STPM.nBins,length(trBatch));
for iter=1:length(trBatch)
    tic;
    thisBlocks = getBlocksPerVid(imdb,trBatch(iter),STPM);
    trFeatures(:,iter) = SCPooling(thisBlocks, dict, STPM.pyramid, SC_opts, LDS_opts);
    t = toc;
    fprintf('Processing training data %d takes %f seconds\n',iter, t);
end
disp('-----------------------------------------------------------End.')

disp('Pooling features for testing data-----------------------------')
teFeatures = zeros(length(dict)*STPM.nBins,length(teBatch));
for iter=1:length(teBatch)
    tic;
    thisBlocks = getBlocksPerVid(imdb,teBatch(iter),STPM);
    teFeatures(:,iter) = SCPooling(thisBlocks, dict, STPM.pyramid, SC_opts, LDS_opts);
    t = toc;
    fprintf('Processing testing data %d takes %f seconds\n',iter, t);
end
disp('-----------------------------------------------------------End.')

disp('Training and predicting')
train_x = trFeatures';
test_x = teFeatures';
train_y = imdb.images.label(trBatch)';
test_y = imdb.images.label(teBatch)';
[accuracy,predict_y] = classify(train_x, train_y, test_x, test_y, Classifier_opts);
disp('-----------------------------------------------------------End.')

t2=clock;
time_dict = etime(t2,t1)

disp('Saving results.')
% save(fullfile(expDir,[num2str(DL_opts.nAtoms),'_STPM_can.mat']),'accuracy', 'dict', 'LDS_opts','SC_opts','DL_opts');