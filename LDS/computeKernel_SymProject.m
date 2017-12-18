function K = computeKernel_SymProject(X1,varargin)
%COMPUTKERNEL: projected kernel under symmetric representation
%  K{i,j} =  inner_product(X1{i},X2{j})

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

n = size(X1{1}.A,1);
K = zeros(numpoints1,numpoints2);


if (same_flag)
    %SY1 = SY2
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2
            
            A1 = X1{t1}.A;
            C1 = X1{t1}.C;
            A2 = X1{t2}.A;
            C2 = X1{t2}.C;
            
            tmpA_t1 = repmat(A1,1,n);
            tmpA_t2 = repmat(A2,1,n);
            E = (1-tmpA_t1.^2).*((1-tmpA_t2'.^2)./(1-tmpA_t1.*tmpA_t2').^2);
            tmpC_12 = C1'*C2;
            tmpKernel_Val = sum(sum((tmpC_12.*E).*tmpC_12));
            
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

            A1 = X1{t1}.A;
            C1 = X1{t1}.C;
            A2 = X2{t2}.A;
            C2 = X2{t2}.C;                       
            
            tmpA_t1 = repmat(A1,1,n);
            tmpA_t2 = repmat(A2,1,n);
            E = (1-tmpA_t1.^2).*((1-tmpA_t2'.^2)./(1-tmpA_t1.*tmpA_t2').^2);
            tmpC_12 = C1'*C2;
            tmpKernel_Val = sum(sum((tmpC_12.*E).*tmpC_12));
            
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
