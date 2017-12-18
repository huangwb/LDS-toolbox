feval('param_default');

%% LDS Parameters
LDS_opts.dim_hidden = 8;
LDS_opts.decomposition = 'CANONICAL';
LDS_opts.Grass_k = 3;
LDS_opts.stabilizer = 'SN';

%% Dictionary Learning Paprameters
DL_opts.nAtoms = 128 ;
DL_opts.vBatchSize = 20;               % 15
DL_opts.nIter = 5;                    % 40
DL_opts.learningRate_A = 0.5;            % 1;0.01
DL_opts.learningRate_C = 5;           % 25;0.5
DL_opts.learningRate_gamma = 0.5;
DL_opts.learningRate_stepsize = 10;

%% Spatial Temporal Pyramid Matching
STPM.blockPerVid = 30;    % 4 16
STPM.bh = 16;             % 8 16 
STPM.bw = 16;             % 8 16
STPM.bt = 50;            % 50
STPM.sh = 16;             % 8 16 
STPM.sw = 16;             % 8 16
STPM.st = 50;            % 50

STPM.pyramid = [1 1 1;
                2 2 2
                4 4 4];

%% Bag-of-System
BoS_opts.computeMean = @computeMean_Grass; 
BoS_opts.nAtoms = 32;
BoS_opts.nIter = 10;
BoS_opts.learningRate_A = 0.5;
BoS_opts.learningRate_C = 1;
BoS_opts.learningRate_gamma = 0.5;
BoS_opts.learningRate_stepsize = 10;

BoS_opts.blockPerVid = 30;
BoS_opts.bh = 16;             % 8 16 
BoS_opts.bw = 16;             % 8 16
BoS_opts.bt = 50;            % 50
BoS_opts.sh = 16;             % 8 16 
BoS_opts.sw = 16;             % 8 16
BoS_opts.st = 50;            % 50

%% Classifier Parameters
% Classifier_opts.classifier = 'LINEARSVM';
Classifier_opts.classifier = 'LIBSVM';
Classifier_opts.kernel = @computeKernel_Chisquare;
feval('param_refresh');


%% Dataset
dbDir = ['../../data/classification/dyntex'];
imgDir = fullfile(dbDir,'gray');
isplit = {'dyntex12'};
expDir = fullfile(pwd, 'output', 'dyntex') ;
setupData = @(x,y,z)setupData_dyntex(x,y,z);
getBlockBatch = @(x,y,z)getBlockBatch_dyntex(x,y,z);
getBlocksPerVid = @(x,y,z)getBlocksPerVid_dyntex(x,y,z);