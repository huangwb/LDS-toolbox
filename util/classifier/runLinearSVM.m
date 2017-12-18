function [accuracy,predict_y] = runLinearSVM(train_x, train_y, test_x, test_y, Classifier_opts)
%RUNNN Summary of this function goes here
%   Detailed explanation goes here

fprintf('Traing classifier using linear SVM.\n');
opts = ['-s ',num2str(1),' -c ',num2str(1)];
classifier.model = train(train_y,sparse(double(train_x)), opts);

fprintf('Predicting.\n');
[predict_y, acc, decValues] = predict(test_y, sparse(double(test_x)), classifier.model);
accuracy = acc(1);

end

