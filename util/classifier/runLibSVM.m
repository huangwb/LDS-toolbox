function [accuracy,predict_y] = runLibSVM(train_x, train_y, test_x, test_y, Classifier_opts)
%RUNNN Summary of this function goes here
%   Detailed explanation goes here

% choose the kernel
computeKernel = Classifier_opts.kernel;

numTrain =  length(train_y);
numTest = length(test_y);

fprintf('Computing the training kernel matrix.\n');
K_train = computeKernel(train_x);
fprintf('Training the svm classifier with the precomputed kernel matrix.\n');
opts = ['-t ' num2str(4) '-b ' num2str(1) '-q ' num2str(1)];
model = ovrtrain(train_y, [(1:numTrain)' K_train], opts);

fprintf('Computing the testing kernel matrix.\n');
K_test = computeKernel(test_x, train_x);

fprintf('Predicting with the pretrained model.\n');
opts = ['-q ' num2str(1)];
[predict_y, predict_accuracy, predict_prob] = ovrpredict(test_y, [(1:numTest)', K_test], model, opts);
accuracy = predict_accuracy*100;

end

