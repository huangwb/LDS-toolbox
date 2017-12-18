function K = computeKernel_Grass(X1,varargin)
%COMPUTKERNEL: projected kernel under Grassmanian representation
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

n = size(X1{1}.O,2);
K = zeros(numpoints1,numpoints2);


if (same_flag)
    %SY1 = SY2
    
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2
            
            
            tmpMatrix = X1{t1}.O'*X1{t2}.O;
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
            
            
            tmpMatrix = X1{t1}.O'*X2{t2}.O;
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
