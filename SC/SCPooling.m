function [features] = SCPooling(blocks, dict, pyramid, SC_opts, LDS_opts)
%================================================
% 
% Usage:
% Compute the linear spatial-temporal pyramid feature using sparse coding. 
%
% Inputss:
% blocks        -structure defining the feature set of an image   
%                   .data      local feature array extracted from the
%                              video, column-wise
%                   .x          x locations of each local feature, 1st
%                               dimension of the matrix
%                   .y          y locations of each local feature, 2nd
%                               dimension of the matrix
%                   .z          z locations of each local feature, 3rd
%                               dimension of the matrix
%                   .width      width of each image frame
%                   .height     height of each image frame
%                   .tlength    number of frames
% dict             -sparse dictionary, column-wise
% pyramid       -defines structure of pyramid 
% 
% Output:
% features          -multiscale max pooling feature
%
% Written by Wenbing Huang

%===============================================
img_width  = blocks.width;
img_height = blocks.height;
img_tlength = blocks.tlength;


nAtoms = size(dict,2);
% nBlocks = features.nBlocks;
% idxBlocks = zeros(nBlocks,1);

% extract lds local feature
lds = getLdsBatch(blocks.data,LDS_opts);

% compute sparse codes for each local feature
sc_codes = sparse_coding(lds,dict,SC_opts);
sc_codes = abs(sc_codes);

% spatial-temperal levels
pLevels = size(pyramid,1);
% spatial bins on each level
pBins = pyramid(:,1).*pyramid(:,2).*pyramid(:,3);
% total spatial bins
tBins = sum(pBins);

features = zeros(nAtoms, tBins);
bId = 0;

for iter1 = 1:pLevels,
    
    nBins = pBins(iter1);
    
    wUnit = img_width / pyramid(iter1,1);
    hUnit = img_height / pyramid(iter1,2);
    tUnit = img_tlength / pyramid(iter1,3); 
    
    % find to which spatial bin each local descriptor belongs
    xBin = ceil(blocks.x / hUnit);
    yBin = ceil(blocks.y / wUnit);
    zBin = ceil(blocks.z / tUnit);
    idxBin = (zBin-1)* pyramid(iter1,2) * pyramid(iter1,2)  + (yBin - 1)*pyramid(iter1,2) + xBin;
    
    for iter2 = 1:nBins,     
        bId = bId + 1;
        sidxBin = find(idxBin == iter2);
        if isempty(sidxBin),
            continue;
        end      
        features(:, bId) = max(sc_codes(:, sidxBin), [], 2);
    end
end

if bId ~= tBins,
    error('Index number error!');
end

features = features(:);
features = features./sqrt(sum(features.^2));
