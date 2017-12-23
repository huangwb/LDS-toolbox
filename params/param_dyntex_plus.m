feval('param_default');

%% LDS Parameters
LDS_opts.dim_hidden = 10;
LDS_opts.decomposition = 'CANONICAL'; % CANONICAL, GRASS, SYMSKEW
LDS_opts.Grass_k = 3;
LDS_opts.stabilizer = 'SN';

%% Dictionary Learning Paprameters
DL_opts.nAtoms = 4;
DL_opts.vBatchSize = 128;             % 15
DL_opts.nIter = 1;                    % 40
DL_opts.learningRate_A = 0.5;            % 1;0.01
DL_opts.learningRate_C = 5;           % 25;0.5
DL_opts.learningRate_gamma = 0.5;
DL_opts.learningRate_stepsize = 10;


%% Clustering Parameters
Clustering_opts.computeMean = @computeMean_Can; 
Clustering_opts.nAtoms = 16;
Clustering_opts.nIter = 20;
Clustering_opts.learningRate_A = 0.1;
Clustering_opts.learningRate_C = 1;
Clustering_opts.learningRate_gamma = 0.5;
Clustering_opts.learningRate_stepsize = 10;

%% Classifier Parameters
Classifier_opts.classifier = 'LINEARSVM';

feval('param_refresh');


%% Dataset
dbDir = ['dataset/dyntex_plus'];
imgDir = fullfile(dbDir,'dyntex++lbp8-6.mat');
isplit = {'1'};
expDir = fullfile(pwd, 'output', 'dyntex++') ;
setupData = @(x,y,z)setupData_memory(x,y,z);
getBlockBatch = @(x,y)getBlockBatch_memory(x,y);
getBlocksPerVid = @(x,y)getBlockBatch_memory(x,y);