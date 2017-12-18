function blocks = getBlocksPerVid_dyntex(imdb,img,sampleOpts)
%GET_BATCH Summary of this function goes here
%   Detailed explanation goes here

% rng('default');
% rng(0);

% load image
imgPath = fullfile(imdb.imageDir, strcat(imdb.images.name{img},'.mat'));
input =  load(imgPath);
video = single((input.frames));

% block statics
blocks.height = size(video,1);
blocks.width = size(video,2);
blocks.tlength = size(video,3);

% hsegs = fix(blocks.height/sampleOpts.bh);
% wsegs = fix(blocks.width/sampleOpts.bw);
% tsegs =  fix(blocks.tlength/sampleOpts.bt);
hsegs = fix((size(video,1)-sampleOpts.bh)/sampleOpts.sh) + 1;
wsegs = fix((size(video,2)-sampleOpts.bw)/sampleOpts.sw) + 1;
tsegs = fix((size(video,3)-sampleOpts.bt)/sampleOpts.st) + 1;

blocks.number = wsegs*hsegs*tsegs;
blocks.data = cell(1,blocks.number);
blocks.x = zeros(1,blocks.number);
blocks.y = zeros(1,blocks.number);
blocks.z = zeros(1,blocks.number);

iter = 1;
for zIter=1:tsegs
    for yIter=1:wsegs
        for xIter=1:hsegs
%             thisBlock = single(video((xIter-1)*sampleOpts.bh+1:xIter*sampleOpts.bh, ...
%                 (yIter-1)*sampleOpts.bw+1:yIter*sampleOpts.bw, ...
%                 (zIter-1)*sampleOpts.bt+1:zIter*sampleOpts.bt));
            thisBlock = single(video((xIter-1)*sampleOpts.sh+1:(xIter-1)*sampleOpts.sh+sampleOpts.bh, ...
                                     (yIter-1)*sampleOpts.sw+1:(yIter-1)*sampleOpts.sw+sampleOpts.bw, ...
                                     (zIter-1)*sampleOpts.st+1:(zIter-1)*sampleOpts.st+sampleOpts.bt));
                                 
            blocks.data{iter} = reshape(thisBlock,sampleOpts.bh*sampleOpts.bw,sampleOpts.bt);
            
%             blocks.x(iter) = (xIter-1)*sampleOpts.bh+1;
%             blocks.y(iter) = (yIter-1)*sampleOpts.bw+1;
%             blocks.z(iter) = (zIter-1)*sampleOpts.bt+1;
            blocks.x(iter) = xIter*sampleOpts.sh+1; 
            blocks.y(iter) = yIter*sampleOpts.sw+1; 
            blocks.z(iter) = zIter*sampleOpts.st+1;
    
            iter = iter + 1;
        end
    end
end

end
