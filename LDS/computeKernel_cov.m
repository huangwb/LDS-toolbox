function K = computeKernel_cov(X1,varargin)
%COMPUTKERNEL: projected kernel between covariance terms 
%   K{i,j} =  inner_product(X1{i},X2{j})

MIN_THRESH = 1e-6;

same_flag = false;
if (nargin < 2)
%     X2 = X1;
    same_flag = true;
    numpoints1 = size(X1,2);
    numpoints2 = numpoints1;
else
    X2 = varargin{:};
    numpoints1 = size(X1,2);
    numpoints2 = size(X2,2);
end

n = size(X1{1}.R,2);
K = zeros(numpoints1,numpoints2);


if (same_flag)
    %SY1 = SY2
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2

            tmpMatrix = (X1{t1}.R)'*(X1{t2}.R);
            
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

             tmpMatrix = (X1{t1}.R)'*(X2{t2}.R);
            
            tmpKernel_Val = sum(sum(tmpMatrix.^2));
            
            if (tmpKernel_Val < MIN_THRESH)
                tmpKernel_Val = 0;
            elseif (tmpKernel_Val > n)
                tmpKernel_Val = n;
            end
%             
            K(t1,t2) = tmpKernel_Val;
        end
    end
end



end

