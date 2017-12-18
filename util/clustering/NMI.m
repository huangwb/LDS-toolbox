function score = NMI(idx_cluster, idx_class)
%PL Summary of this function goes here
%   Detailed explanation goes here

nData = length(idx_cluster);
nClusters = length(unique(idx_cluster));
nClasses =  length(unique(idx_class));

count_class = hist(idx_class,unique(idx_class));
count_cluster = hist(idx_cluster,unique(idx_cluster));


% compute mutual information
I = 0;
for k=1:nClusters
    this_cluster = (idx_cluster==k);
    count_class_k = hist(idx_class(this_cluster),1:nClasses);
    valid_class_k = (count_class_k>0);
    I = I + sum(1/nData*count_class_k(valid_class_k).*log(nData/count_cluster(k)*count_class_k(valid_class_k)./count_class(valid_class_k)));
end

% compute the entropy for clustering
H_cluster = -sum(1/nData*count_cluster.*log(1/nData*count_cluster));
% compute the entropy for classification
H_class = -sum(1/nData*count_class.*log(1/nData*count_class));


score = 2*I/(H_cluster+H_class);

end

