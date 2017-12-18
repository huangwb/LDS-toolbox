function K = computeKernel_Martin(X1,varargin)
%COMPUTKERNEL: RBF kernel with Martin distance
%   K{i,j} =  inner_product(X1{i},X2{j})

p = inputParser;
sigma2 =100;
p.addOptional('parameter',sigma2);
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
    case 3
        same_flag = true;
        numpoints1 = size(X1,2);
        numpoints2 = numpoints1;
        parse(p,varargin{:});
        sigma2 = p.Results.parameter;
    case 4
        same_flag = false;
        X2 = varargin{1};
        numpoints1 = size(X1,2);
        numpoints2 = size(X2,2);
        parse(p,varargin{:});
        sigma2 = p.Results.parameter;
end

K = zeros(numpoints1,numpoints2);

if (same_flag)
    %SY1 = SY2
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2
            
            d = distMartin(X1{t1}.A,X1{t1}.C,X1{t2}.A,X1{t2}.C);
            tmpKernel_Val = exp(-d/sigma2);
            
            
            K(t1,t2) = tmpKernel_Val;
            K(t2,t1) = K(t1,t2);
            
        end
    end
else
    for t1 = 1:numpoints1
        for t2 = 1:numpoints2
            
            d = distMartin(X1{t1}.A,X1{t1}.C,X2{t2}.A,X2{t2}.C);
            tmpKernel_Val = exp(-d/sigma2);
            
            K(t1,t2) = tmpKernel_Val;
        end
    end
end

end