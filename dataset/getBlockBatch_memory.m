function blocks = getBlockBatch_memory(imdb,imBatch,sampleOpts)
%GET_BATCH gets blockbatch from memory
%   Detailed explanation goes here
blocks.data = imdb.images.frames(imBatch);

end

