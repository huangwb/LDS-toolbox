function [accuracy,predict_y] = runNN(train_x, train_y, test_x, test_y, Classifier_opts)
%RUNNN Summary of this function goes here
%   Detailed explanation goes here

% choose the kernel 
computeKernel = Classifier_opts.kernel; 

fprintf('Computing the kernel matrix.\n');
D = computeKernel(train_x,test_x);
[~,maxIndex] = max(D);
predict_y = train_y(maxIndex);

accuracy = sum(predict_y == test_y)/length(test_y)*100;

end

