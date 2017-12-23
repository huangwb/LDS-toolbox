%% LDS Parameters
% basic parameters 
LDS_opts.dim_hidden = 10;                       % The dimensionality of the hidden state (i.e. n).            
LDS_opts.nv = 1;                                % The dimensionality of the covariance matrix.
LDS_opts.transition = 'FREE';                   % Two kinds of transition matrices can be selected: SYM or FREE
LDS_opts.decomposition = 'NONE';                % Four choice: NONE; CANONICAL; SYMSKEW; GRASS
LDS_opts.Grass_k = 5;                           % The time length for the Grass decomposition
% stabilization parameters.
LDS_opts.stabilizer = 'NONE';         % Methods for stabilization. none: no stabilization; SN: Soft normalization with sigmoid function; UN: unified-normalization; MN: maximized-normalization; HZP: hard-zero-padding; WLS: weighted-least-square; DWLS: digonal-weighted-least-square
                                      % CG: constraint-geneartion; LB1; LB2;
LDS_opts.SN_a = 2.5;                  % Parameters for SN
LDS_opts.SN_b = 1;
LDS_opts.UN_factor = 1.5;             % Only valid for the UN method.
LDS_opts.bound = 1;                   % The bound value of the Spectral Radius of the transition matrix, only valid for WLS, DWLS, CG, LB2
LDS_opts.svd_check = 0;               % Setting 0 means checking the bounds of eigent values, while 1 means checking the bounds of singular values. Only valid for WLS, DWLS, CG

LDS_opts.kernel = @computeKernel_Project;       % See the definitions of the following kernels in ./LDS/ : @computeKernel_Project; @computeKernel_CanProject; @computeKernel_SymProject; @computeKernel_SymSkew; @computeKernel_Martin
LDS_opts.beta = 1;                              % The weight determining the trade-off between the mean and covariance components. Setting 1 means no covariance is involved, while 0 means no mean is involved. 


%% Sparse Coding Parameters
SC_opts.lambda = 0.1;                          % Sparsity penalty weight
SC_opts.kernel = LDS_opts.kernel;
SC_opts.L = 10;                                % The number of dictionary atoms for local sparse coding


%% Dictionary Learning Paprameters
DL_opts.nAtoms = 15;
DL_opts.vBatchSize = 15;               % 15
DL_opts.nIter = 40;                    % 40
DL_opts.learningRate_A = 1;            % 1;0.01
DL_opts.learningRate_C = 25;           % 25;0.5
DL_opts.learningRate_gamma = 0.5;
DL_opts.learningRate_stepsize = 50;
DL_opts.initial_method = 'RAND';
DL_opts.teBatchSize = 40;


%% Clustering Parameters
Clustering_opts.kernel = LDS_opts.kernel;
Clustering_opts.computeMean = @computeMean_MDS; % methods applied for mean computation of LDS data: @computeMean_MDS and @computeMean_Can
Clustering_opts.nAtoms = 15;
Clustering_opts.nIter = 5;
Clustering_opts.learningRate_A = 0.1;
Clustering_opts.learningRate_C = 1;


%% Classifier Parameters
Classifier_opts.classifier = 'SCC';  % methods applied for classification: 'SCC', 'NN', 'SVM', 'LINEARSVM'
Classifier_opts.kernel = LDS_opts.kernel;
Classifier_opts.SC_opts = SC_opts;

