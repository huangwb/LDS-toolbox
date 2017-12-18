function K = computeKernel_CanProject(X1,varargin)
%COMPUTKERNEL: projected kernel under canonical representation
%   K{i,j} =  inner_product(X1{i},X2{j})

MIN_THRESH = 1e-6;

switch nargin
    case 1
        same_flag = true;
        numpoints1 = size(X1,2);
        numpoints2 = numpoints1;
    case 2
        same_flag = false;
        X2 = varargin{1};
        numpoints1 = size(X1,2);
        numpoints2 = size(X2,2);
end

n = size(X1{1}.S,1);
K = zeros(numpoints1,numpoints2);


if (same_flag)
    %SY1 = SY2
    
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2
            
            A1 = diag(X1{t1}.S)*X1{t1}.V;
            C1 = X1{t1}.U;
            A2 = diag(X1{t2}.S)*X1{t2}.V;
            C2 = X1{t2}.U;
            
            iL1 = X1{t1}.iL;
            iL2 = X1{t2}.iL;
%             iL1 = compute_invL(A1);
%             iL2 = compute_invL(A2);
            
            tmpMatrix = iL1 * dlyap(A1', A2, C1'*C2)* iL2';
            tmpKernel_Val = sum(sum(tmpMatrix.^2));
            
            if (tmpKernel_Val < MIN_THRESH)
                tmpKernel_Val = 0;
            elseif (tmpKernel_Val > n)
                tmpKernel_Val = n;
            end
            
            K(t1,t2) = tmpKernel_Val;
            K(t2,t1) = K(t1,t2);
        end
    end
else
    for t1 = 1:numpoints1
        for t2 = 1:numpoints2
            
            A1 = diag(X1{t1}.S)*X1{t1}.V;
            C1 = X1{t1}.U;
            A2 = diag(X2{t2}.S)*X2{t2}.V;
            C2 = X2{t2}.U;     
            
            iL1 = X1{t1}.iL;
            iL2 = X2{t2}.iL;
%             iL1 = compute_invL(A1);
%             iL2 = compute_invL(A2);
            
            tmpMatrix = iL1 * dlyap(A1', A2, C1'*C2)* iL2';
            tmpKernel_Val = sum(sum(tmpMatrix.^2));
            
            if (tmpKernel_Val < MIN_THRESH)
                tmpKernel_Val = 0;
            elseif (tmpKernel_Val > n)
                tmpKernel_Val = n;
            end
            
            K(t1,t2) = tmpKernel_Val;
        end
    end
end

end

% function invL = compute_invL(A)
% n = size(A,1);
% dot = dlyap(A',eye(n));
% [Udot,Sdot,~] = svd(dot);
% invL = diag(1./sqrt(diag(Sdot)))*Udot';
% end