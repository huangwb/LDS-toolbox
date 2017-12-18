function [s,P] = symTranMatrix(X0,X1)
% symTranMatrix returns symmetric transition matrix
%
% INPUTS
% X0         - [] Current state vetors
% X1         - [] One-time-step forward state vectors
%
% OUTPUTS
% s          - [] Diagonal elements of the standar form of the sysmmetric transition matrix
% P          - [] Othorgonal matrix of the standar form of the sysmmetric transition matrix
%
% implemented by Wenbing Huanng, 2016-4-04

n = size(X0,1);

A = X1*pinv(X0);
[U,S,~]=svd(1/2*(A+A'));
tmp = diag(S);
s = tmp(1:n);
P = U(:,1:n);


nIters = 5;

for k=1:n-1
    if k>1
        W = null(P(:,1:k-1)');
    else
        W = eye(n);
    end   
    
    for iter=1:nIters
        s(k) = ((P(:,k)'*X0)*(X1'*P(:,k)))/((P(:,k)'*X0)*(X0'*P(:,k)));
        M = s(k)^2*(W'*X0)*(X0'*W)- s(k)*(W'*X0)*(X1'*W)-s(k)*(W'*X1)*(X0'*W);
        eigVector=smallestEigvector(M,1);
%         eigVector = GCG( W'*P(:,k),M,20 );
        P(:,k) = W*eigVector;
    end
end

P(:,n) =  null(P(:,1:n-1)');
s(n) = ((P(:,n)'*X0)*(X1'*P(:,n)))/((P(:,n)'*X0)*(X0'*P(:,n)));

end

