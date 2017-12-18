function [eigvector,eigenvalue]= smallestEigvector(S,k)
% smallestEigvector returns the smallest eigenvectors and eigenvalues with
% MATLAB embedded function
%
% INPUTS
% S         - [n n] the scale matrix of the SVD for observation vectors
% k         - [1] number of the eigenvectors needed to return
%
% OUTPUTS
% eigvector      - [n,k] the resulting eigenvectors
% eigenvalue     - [k,k] the resulting eigenvalues
%
% implemented by Wenbing Huanng, 2016-4-04
S = double(real(S));
S = 1/2*(S+S');
% [v,maxlambda] = eigs(S,1);
% if maxlambda<0
%     eigvector = v;
% else
%     [eigvector,maxlambda] = eigs(maxlambda*eye(size(S))-S,1);
% end
opts.issym = 1;
opts.isreal = 1;
[eigvector,~]=eigs(S,k,'SA',opts);
eigvector = real(eigvector);
eigenvalue = eigvector'*S*eigvector;

end