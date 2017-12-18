function score = purity(idx_cluster, idx_class)
%PL Summary of this function goes here
%   Detailed explanation goes here

nData = length(idx_cluster);
nClusters = length(unique(idx_cluster));
nClasses =  length(unique(idx_class));

score = 0;
for k=1:nClusters
    this_cluster = (idx_cluster==k);
    count_class_k = hist(idx_class(this_cluster),1:nClasses);
    score = score + max(count_class_k);
end

score = score/nData;

end

