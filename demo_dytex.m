% Here is an example to learn the optimal LDS tuples for the input video.
% The default learner is the subspace method developed by Doretto et al.
% To apply other learning methods, you need to set the parameter LDS_opts. 
% See ./params/default.m for more details.

clear
close
addPaths; % Add all ducuments to the current path

disp('Set-ups')
feval('param_dyntex_plus');

disp('Preparing dataset')
imdbPath = fullfile(expDir, sprintf('imdb_%s.mat',isplit{1}));
if exist(imdbPath, 'file')
    imdb = load(imdbPath) ;
else
    imdb = setupData(dbDir,imgDir,isplit{1});
    if ~exist(expDir, 'dir')
        mkdir(expDir) ;
    end
    save(imdbPath, '-struct', 'imdb') ;
end

thisBlocks = getBlocksPerVid(imdb,1);
thisLds = getLdsBatch(thisBlocks.data,LDS_opts); 

% if you want to lear a stable system, you can apply the weighted square
% method proposed by our IJCAI_2016 paper
LDS_opts.stabilizer = 'WLS';
thisLds = getLdsBatch(thisBlocks.data,LDS_opts); 

