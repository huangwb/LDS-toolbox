function [Q, d] = align_Group(X1,X2,opts)
%ALIGN_GROUP Summary of this function goes here
%   Detailed explanation goes here

% params
lambda_C =  1;
A1 = diag(X1.S)*X1.V;
C1 = X1.U;
A2 = diag(X2.S)*X2.V;
C2 = X2.U;
n = size(A1,1);

% random initialization
pos_Q = orth(randn(n));
% in case rank(pos_Q)<n
res = null(pos_Q');
if ~isempty(res)
    pos_Q = [pos_Q;res];
end

neg_Q = -pos_Q;

% align X1 to X2
for iter=1:opts.nIter
    % compute the gradient w.r.t. A
    delta_pos_Q = A1'*A1*pos_Q - A1'*pos_Q*A2 + A1*A1'*pos_Q - A1*pos_Q*A2';
    delta_neg_Q = A1'*A1*neg_Q - A1'*neg_Q*A2 + A1*A1'*neg_Q - A1*neg_Q*A2';
    
    % compute the gradient w.r.t. C
    delta_pos_Q = delta_pos_Q + lambda_C * (pos_Q - C1'*C2);
    delta_neg_Q = delta_neg_Q + lambda_C * (neg_Q - C1'*C2);
    
    % update via Caley-transformation
    pos_Q = project_to_manifold(pos_Q, delta_pos_Q, opts.learningRate_A);
    neg_Q = project_to_manifold(neg_Q, delta_neg_Q, opts.learningRate_A);
    
end

f_pos_Q = sum(sum((pos_Q'*A1*pos_Q - A2).^2)) + lambda_C * sum(sum((C1*pos_Q - C2).^2  ));
f_neg_Q = sum(sum((neg_Q'*A1*neg_Q - A2).^2)) + lambda_C * sum(sum((C1*neg_Q - C2).^2  ));

Q = pos_Q;
d = f_pos_Q;
if f_neg_Q<f_pos_Q
    Q = neg_Q;
    d = f_neg_Q;
end


end

