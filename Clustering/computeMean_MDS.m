function center = computeMean_MDS(data, opts, varargin)
%COMPUTMEAN_MDS Summary of this function goes here
%   Detailed explanation goes here

n = length(data{1}.S);
nData = length(data);
computeKernel = opts.kernel;

% compute kernel matrix
% distM = 2*(n - computeKernel(data))/n;
% distM = distM - diag(diag(distM));
distM = computeKernel(data)/n;
distM = distM - diag(diag(distM)) + diag(ones(nData,1));

% perform MultiDimensional Scaling (MDS)
Y = cmdscale(distM);

% compute the mean in the projected space
Y_mean = mean(Y);

% find the nearest data to the mean in the projected space
dist2mean = sum(bsxfun(@minus,Y,Y_mean).^2,2);
[~, id] = min(dist2mean); 

center = data{id};

end

