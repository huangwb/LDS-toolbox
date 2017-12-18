function K = computeKernel_Chisquare(X1,varargin)
%COMPUTKERNEL: RBF kernel with Chi-Square distance
%   K{i,j} =  inner_product(X1(i,:),X2(j,:))

p = inputParser;
sigma2 =1;
p.addOptional('parameter',sigma2);
switch nargin
    case 1
        same_flag = true;
        numpoints1 = size(X1,1);
        numpoints2 = numpoints1;
    case 2
        same_flag = false;
        X2 = varargin{1};
        numpoints1 = size(X1,1);
        numpoints2 = size(X2,1);
    case 3
        same_flag = true;
        numpoints1 = size(X1,1);
        numpoints2 = numpoints1;
        parse(p,varargin{:});
        sigma2 = p.Results.parameter;
    case 4
        same_flag = false;
        X2 = varargin{1};
        numpoints1 = size(X1,1);
        numpoints2 = size(X2,1);
        parse(p,varargin{:});
        sigma2 = p.Results.parameter;
end

K = zeros(numpoints1,numpoints2);

if (same_flag)
    %SY1 = SY2
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2
            valid = ((X1(t1,:)+X1(t2,:))>0);
            d = 1/2*sum((X1(t1,valid)-X1(t2,valid)).^2./(X1(t1,valid)+X1(t2,valid)));
            tmpKernel_Val = exp(-d/sigma2);
            
            
            K(t1,t2) = tmpKernel_Val;
            K(t2,t1) = K(t1,t2);
            
        end
    end
else
    for t1 = 1:numpoints1
        for t2 = 1:numpoints2
            valid = ((X1(t1,:)+X2(t2,:))>0);
            d = 1/2*sum((X1(t1,valid)-X2(t2,valid)).^2./(X1(t1,valid)+X2(t2,valid)));
            tmpKernel_Val = exp(-d/sigma2);
            
            K(t1,t2) = tmpKernel_Val;
        end
    end
end

end