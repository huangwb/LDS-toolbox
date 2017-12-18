function [accuracy, predict_y] = SCC( train_x, train_y, test_x, test_y, Classifier_opts)
%SCC classification with reconstruction error approach
%   Detailed explanation goes here

fprintf('Computing codes for testing data.\n');
[alpha,D,qX]=sparse_coding(test_x,train_x,Classifier_opts.SC_opts);

fprintf('Predicting.\n');
numClass = max(train_y);
res_y = zeros(numClass,size(alpha,2));
for t1 = 1:numClass
    classIndex = (train_y==t1);
    delta_alpha = zeros(size(alpha));
    delta_alpha(classIndex,:) = alpha(classIndex,:);
    res_y(t1,:) = sum((qX - D*delta_alpha).^2);
end

[~,minIndex] = min(res_y);
predict_y = minIndex';

accuracy = sum(predict_y == test_y)/length(test_y)*100;


end

