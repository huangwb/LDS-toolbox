function K = computeKernel_SymSkew(X1,varargin)
%COMPUTKERNEL: projected kernel under two-fold representation
%    K{i,j} =  inner_product(X1{i},X2{j})

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

n = size(X1{1}.A1,1);
K = zeros(numpoints1,numpoints2);


% if isfield(X1{1}, 'A')
%     n = size(X1{1}.A,1);
% elseif isfield(X1{1}, 'A1')
%     n = size(X1{1}.A1,1);
% end

compute_kernel_function = @compute_kernel_val_22;

if (same_flag)
    %SY1 = SY2
    for t1 = 1:numpoints1
        for t2 = t1:numpoints2            
            
            kernel_val = compute_kernel_function(X1{t1}, X1{t2}, n);
            
            if (kernel_val < MIN_THRESH)
                kernel_val = 0;
            elseif (kernel_val > n)
                kernel_val = n;
            end
            K(t1,t2) = kernel_val;
            K(t2,t1) = K(t1,t2);
        end
    end
else
    for t1 = 1:numpoints1
        for t2 = 1:numpoints2
            
            kernel_val = compute_kernel_function(X1{t1}, X2{t2}, n);
            
            if (kernel_val < MIN_THRESH)
                kernel_val = 0;
            elseif (kernel_val > n)
                kernel_val = n;
            end
            K(t1,t2) = kernel_val;
        end
    end
end

end

% function kernel_val = compute_kernel_val_11(X1, X2,n)
% 
% A1 = X1.A;
% C1 = X1.C;
% A2 = X2.A;
% C2 = X2.C;
% 
% tmpD1 = repmat(A1,1,n);
% tmpD2 = repmat(A2,1,n);
% 
% E12 = (1-tmpD1.*conj(tmpD1)).*(1-tmpD2.*conj(tmpD2))'./((1-tmpD1.*tmpD2').*conj(1-tmpD1.*tmpD2'));
% tmpC12 = (C1)'*(C2);
% 
% kernel_val = abs(sum(sum((tmpC12.*E12).*conj(tmpC12))));
% 
% end
% 
% function kernel_val = compute_kernel_val_12(X1, X2,n)
% 
% A1 = X1.A;
% C1 = X1.C;
% A2_1 = X2.A1;
% A2_2 = X2.A2;
% C2_1 = X2.C1;
% C2_2 = X2.C2;
% 
% tmpD1   = repmat(A1,1,n);
% tmpD2_1 = repmat(A2_1,1,n);
% tmpD2_2 = repmat(A2_2,1,n);
% 
% E12_11 = (1-tmpD1.*conj(tmpD1)).*(1-tmpD2_1.*conj(tmpD2_1))'./((1-tmpD1.*tmpD2_1').*conj(1-tmpD1.*tmpD2_1'));
% E12_22 = (1-tmpD1.*conj(tmpD1)).*(1-tmpD2_2.*conj(tmpD2_2))'./((1-tmpD1.*tmpD2_2').*conj(1-tmpD1.*tmpD2_2'));
% 
% tmpC12_11 = (C1)'*(C2_1);
% tmpC12_22 = (C1)'*(C2_2);
% 
% kernel_val = 1/2*abs(sum(sum((tmpC12_11.*E12_11).*conj(tmpC12_11))) + sum(sum((tmpC12_22.*E12_22).*conj(tmpC12_22))));
% 
% end
% 
% function kernel_val = compute_kernel_val_21(X1,X2,n)
% 
% kernel_val = compute_kernel_val_21(X2,X1,n);
% 
% end

function kernel_val = compute_kernel_val_22(X1, X2,n)

A1_1 = X1.A1;
A1_2 = X1.A2;
C1_1 = X1.C1;
C1_2 = X1.C2;
A2_1 = X2.A1;
A2_2 = X2.A2;
C2_1 = X2.C1;
C2_2 = X2.C2;

tmpD1_1 = repmat(A1_1,1,n);
tmpD1_2 = repmat(A1_2,1,n);
tmpD2_1 = repmat(A2_1,1,n);
tmpD2_2 = repmat(A2_2,1,n);

E12_11 = (1-tmpD1_1.*conj(tmpD1_1)).*(1-tmpD2_1.*conj(tmpD2_1))'./((1-tmpD1_1.*tmpD2_1').*conj(1-tmpD1_1.*tmpD2_1'));
E12_22 = (1-tmpD1_2.*conj(tmpD1_2)).*(1-tmpD2_2.*conj(tmpD2_2))'./((1-tmpD1_2.*tmpD2_2').*conj(1-tmpD1_2.*tmpD2_2'));

tmpC12_11 = (C1_1)'*(C2_1);
tmpC12_22 = (C1_2)'*(C2_2);

kernel_val = 1/2*abs(sum(sum((tmpC12_11.*E12_11).*conj(tmpC12_11))) + sum(sum((tmpC12_22.*E12_22).*conj(tmpC12_22))));

end