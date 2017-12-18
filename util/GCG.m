function [Y,F] = GCG( Y0,A,iters )
%GCG Grassmaniann-Conjungate-Gradient method to find the smallest
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

n = size(Y0,1);
Y = Y0;
G =  A*Y - Y*(Y'*A*Y);
H = -G;
for k=1:iters
    s = sqrt(H'*H);
    U = 1/s*H;
    t = 0.1;
    
    Y1 = Y*cos(t)+U*sin(t);
    G1 = A*Y1-Y1*(Y1'*A*Y1);
    
    tH = (-Y*sin(t)+U*cos(t))*s;
    tG = G - (Y*sin(t)+U*(1-cos(t)))*(U'*G);
    
    gama = ((G1-tG)'*G1)/(G'*G);
    H = -G1 + gama*tH;
    
    if mod(k+1,n-1)
        H = -G1;
    end
    Y = Y1;
end
Y = 1/(Y'*Y)*Y;
F = Y'*A*Y;
% t = toc
end

