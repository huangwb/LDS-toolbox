function codes = BoSCoding(data, center, opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

computeKernel = opts.kernel;

l_sim = computeKernel(center,data);
[~,cluster_idx] = max(l_sim);

codes = hist(cluster_idx,1:opts.nAtoms)/length(cluster_idx);

end

