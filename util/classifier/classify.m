function [accuracy,predict_y] = classify(train_x, train_y, test_x, test_y, Classifier_opts)
%CLASSIFY performs classification 
%   Detailed explanation goes here

switch Classifier_opts.classifier
    case 'NN'
        [accuracy,predict_y] = runNN(train_x, train_y, test_x, test_y, Classifier_opts);
    case 'SCC'
        [accuracy,predict_y] = SCC(train_x, train_y, test_x, test_y, Classifier_opts);
    case 'LINEARSVM'
        [accuracy,predict_y] = runLinearSVM(train_x, train_y, test_x, test_y, Classifier_opts);
    case 'LIBSVM'
        [accuracy,predict_y] = runLibSVM(train_x, train_y, test_x, test_y, Classifier_opts);
end

end

