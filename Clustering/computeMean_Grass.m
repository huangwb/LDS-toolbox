function center = computeMean_Grass(data, opts, varargin)
%COMPUTMEAN_MDS Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(data{1}.O);
nData = length(data);

tmpBig = zeros(m,m);
for i = 1:nData
    tmpBig = tmpBig + data{i}.O*data{i}.O';
end

[O,~] = eigs(1/2*double(tmpBig+tmpBig'),n);
center.O = O;

end

