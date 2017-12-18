function [accuracy, expect_labels] = ridgeRegression(train_alpha, train_y, test_alpha, test_y)
%Ridge regression

num_classes = max(train_y);
num_data = length(train_y);
L = zeros(num_classes,num_data);
L(sub2ind([num_classes,num_data], train_y', 1:num_data)) = 1;
zeta = 1e-1;    %regularization parameter for ridge regression
expect_alpha = train_alpha*train_alpha' + zeta*eye(size(train_alpha,1));
expect_v = L*train_alpha';
W = expect_v/expect_alpha;

[~,expect_labels] = max(W*test_alpha);  

accuracy = sum(expect_labels == test_y')/length(expect_labels);

end