function blocks = getBlockBatch_dyntex(imdb,imBatch,sampleOpts)
%GET_BATCH Summary of this function goes here
%   Detailed explanation goes here

% rng('default');
% rng(0);

numImBatch = numel(imBatch);
blocks.number = numImBatch * sampleOpts.blockPerVid;

blocks.data = cell(1,blocks.number);
blocks.x = zeros(1,blocks.number);
blocks.y = zeros(1,blocks.number);
blocks.z = zeros(1,blocks.number);

for i=1:numImBatch
    
    input =  load(fullfile(imdb.imageDir, strcat(imdb.images.name{imBatch(i)},'.mat')));
    video = single((input.frames));
   
%     hsegs = fix(size(video,1)/sampleOpts.bh);
%     wsegs = fix(size(video,2)/sampleOpts.bw);
%     tsegs =  fix(size(video,3)/sampleOpts.bt);
    hsegs = fix((size(video,1)-sampleOpts.bh)/sampleOpts.sh) + 1;
    wsegs = fix((size(video,2)-sampleOpts.bw)/sampleOpts.sw) + 1;  
    tsegs = fix((size(video,3)-sampleOpts.bt)/sampleOpts.st) + 1;  
      
      
    if wsegs*hsegs*tsegs<sampleOpts.blockPerVid
        error('The number of sampling blocks should be less than that of all possible block segmentations for each image.');
    end    
    
    x = randi(hsegs,[1,sampleOpts.blockPerVid]);
    y = randi(wsegs,[1,sampleOpts.blockPerVid]);
    z = randi(tsegs,[1,sampleOpts.blockPerVid]);
    
    for j=1:sampleOpts.blockPerVid
        
%         thisBlock = single(video((x(j)-1)*sampleOpts.bh+1:x(j)*sampleOpts.bh, ...
%                                  (y(j)-1)*sampleOpts.bw+1:y(j)*sampleOpts.bw, ...
%                                  (z(j)-1)*sampleOpts.bt+1:z(j)*sampleOpts.bt));

%         if (x(j)-1)*sampleOpts.sh+sampleOpts.bh>size(video,1) ||  (y(j)-1)*sampleOpts.sw+sampleOpts.bw > size(video,2) || (z(j)-1)*sampleOpts.st+sampleOpts.bt > size(video,3)
%             stop = 1;
%         end
        
        thisBlock = single(video((x(j)-1)*sampleOpts.sh+1:(x(j)-1)*sampleOpts.sh+sampleOpts.bh, ...
                                 (y(j)-1)*sampleOpts.sw+1:(y(j)-1)*sampleOpts.sw+sampleOpts.bw, ...
                                 (z(j)-1)*sampleOpts.st+1:(z(j)-1)*sampleOpts.st+sampleOpts.bt));                             
                             
                             
        blocks.data{(i-1)*sampleOpts.blockPerVid+j} = reshape(thisBlock,sampleOpts.bh*sampleOpts.bw,sampleOpts.bt);
        
    end
    
%     blocks.x((i-1)*sampleOpts.blockPerVid+1:i*sampleOpts.blockPerVid) = x*sampleOpts.bh+1; 
%     blocks.y((i-1)*sampleOpts.blockPerVid+1:i*sampleOpts.blockPerVid) = y*sampleOpts.bw+1; 
%     blocks.z((i-1)*sampleOpts.blockPerVid+1:i*sampleOpts.blockPerVid) = z*sampleOpts.bt+1;
    blocks.x((i-1)*sampleOpts.blockPerVid+1:i*sampleOpts.blockPerVid) = x*sampleOpts.sh+1; 
    blocks.y((i-1)*sampleOpts.blockPerVid+1:i*sampleOpts.blockPerVid) = y*sampleOpts.sw+1; 
    blocks.z((i-1)*sampleOpts.blockPerVid+1:i*sampleOpts.blockPerVid) = z*sampleOpts.st+1;    
    
end

% save the dimensionalities of the final video for future application
blocks.height = size(video,1);
blocks.width = size(video,2);
blocks.tlength = size(video,3);

end

