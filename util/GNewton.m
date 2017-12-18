function [Y,F] = GNewton(Y,A,iters)
%GNewton Grassmaniann-Newton method to find the smallest
%eigenvalue and eigenvector
%
% INUTS
% Y0           - [] Initial eigenvector 
% A            - [] Target matrix
% iters        - [] Number of iterations
%
% OUTPUTS
% Y          - [] The resulted eigenvector
% F          - [] The resulted eigenvalue
%
% implemented by Wenbing Huanng, 2016-4-04
% tic


n = size(Y,1);
Y0 = Y;
for i=1:5
    delta = -A*Y0 + Y0*(Y0'*A*Y0);
    s = sqrt(delta'*delta);
    U = 1/s*delta;
    s = 0.1;
    Y = Y0*cos(s)+U*sin(s);
    Y0 = Y;    
end

for i=1:iters
    delta = (A-Y0*(Y0'*A)-Y0'*A*Y0*eye(n))\((Y0'*A*Y0)*Y0-A*Y0);
    s = sqrt(delta'*delta);
    U = 1/s*delta;
    Y = Y0*cos(s)+U*sin(s);
    Y0 = Y;
end
F = Y'*A*Y;
end

