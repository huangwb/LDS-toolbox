function score = Fmeasure(idx_cluster, idx_class)
%PL Summary of this function goes here
%   Detailed explanation goes here

nData = length(idx_cluster);
nClusters = length(unique(idx_cluster));
nClasses =  length(unique(idx_class));

count_class = hist(idx_class,1:nClasses);
count_cluster = hist(idx_cluster,1:nClusters);





end

