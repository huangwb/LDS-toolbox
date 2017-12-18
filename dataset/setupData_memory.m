function imdb = setupData_memory(dbDir,imgDir,isplit)
%SETUPDATA Summary of this function goes here
%   Detailed explanation goes here

rng('default');
rng(0);
input = load(imgDir);
trainind = input.trainind;
testind = input.testind;
imdb.imageDir = imgDir;
imdb.images.label = input.label' ;
imdb.images.frames = cell(1,numel(input.label));
for iter=1:numel(input.label)
%     imdb.images.frames{iter} = reshape(input.data{iter},48,48,size(input.data{iter},2));
    imdb.images.frames{iter} = input.data{iter};
end
imdb.images.id = 1:numel(input.label);
imdb.images.set = zeros(1, numel(input.label)) ;
imdb.images.set(trainind{str2double(isplit)}) = 1;
imdb.images.set(testind{str2double(isplit)}) = 3;

end

