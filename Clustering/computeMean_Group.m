function [center,d] = computeMean_Group(data, opts, varargin)
%COMPUTMEAN_MDS Summary of this function goes here
%   Detailed explanation goes here

n = length(data{1}.S);
nData = length(data);
computeKernel = opts.kernel;

% compute kernel matrix
distM = computeKernel(data)/n;
distM = distM - diag(diag(distM)) + diag(ones(nData,1));

% initialization
[~,minIdx] = max(sum(distM));
% minIdx = randi(nData,1);
center = data{minIdx};
Q = cell(1,nData);
d = zeros(1,nData);
for iter=1:opts.nIter
    
    % align data to center
    sum_A = 0;
    sum_C = 0;
    for id_data=1:nData
        
        [this_Q, this_d] = align_Group(data{id_data},center,opts);
        Q{id_data} = this_Q;
        d(id_data) = this_d;
        
        sum_A = sum_A + this_Q'* (diag(data{id_data}.S)*data{id_data}.V)*this_Q;
        sum_C = sum_C + data{id_data}.U * this_Q;
    end
    
    sum_A = 1/nData*sum_A;
    sum_C = 1/nData*sum_C;
    
    % update center
    [U_A,S_A,V_A] = svd(sum_A);
    [U_C,S_C,V_C] = svd(sum_C,'econ');
    center.S = diag(S_A);
    center.V = V_A'*U_A;
    center.U = U_C*V_C'*U_A;
    % update iL
    A_r = diag(center.S)*center.V; 
    S_rr = dlyap(A_r',eye(n));
    [Udot,Sdot,~] = svd(S_rr);
    center.iL = diag(1./sqrt(diag(Sdot)))*Udot';
end

end

