clear
close
addPaths; % Add all ducuments to the current path

rng('default');
rng(0);

disp('Set-ups')
feval('param_dyntex_plus');

disp('Preparing dataset')
imdbPath = fullfile(expDir, sprintf('imdb_%s.mat',num2str(isplit{1})));
if exist(imdbPath, 'file')
    imdb = load(imdbPath) ;
else
    imdb = setupData(dbDir,imgDir,isplit{1});
    if ~exist(expDir, 'dir')
        mkdir(expDir) ;
    end
    save(imdbPath, '-struct', 'imdb') ;
end
trBatch = find(imdb.images.set==1);
teBatch = find(imdb.images.set==3);
batch = 1:length(imdb.images.set);

score_NMI = zeros(1,length(BoS_opts.nAtoms));
score_purity = zeros(1,length(BoS_opts.nAtoms));
time = zeros(1,length(BoS_opts.nAtoms));
accuracy = zeros(1,length(BoS_opts.nAtoms));

for t=1:length(BoS_opts.nAtoms)
    tBoS_opts = BoS_opts;
    tBoS_opts.nAtoms =BoS_opts.nAtoms(t);
    
    t1=clock;
    disp('Clustering via K-means-----------------------------')
    blocks = getBlockBatch(imdb,batch,tBoS_opts);
    lds = getLdsBatch(blocks.data,LDS_opts);
    [centers,idx_cluster] = Kmeans(lds, tBoS_opts);
    disp('-----------------------------------------------End.')
    t2=clock;
    time(t) = etime(t2,t1);
    
    disp('Evaluating clustering scores')
    score_NMI(t) = NMI(idx_cluster, imdb.images.label(batch));
    score_purity(t) = purity(idx_cluster, imdb.images.label(batch));
    disp('-----------------------------------------------------------End.')
    
    fprintf('NMI: %f; Purity: %f; Time: %f; nAtoms: %d\n',score_NMI(t), score_purity(t), time(t), BoS_opts.nAtoms(t));
    
end
% save results
% disp('Saving results.')
% save(fullfile(expDir,['Clustering_Can.mat']),'score_NMI', 'score_purity', 'time');
