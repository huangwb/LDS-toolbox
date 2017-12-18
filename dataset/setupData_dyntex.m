function imdb = setupData_dyntex(dbDir,imgDir,isplit)
%SETUPDATA Summary of this function goes here
%   Detailed explanation goes here
rng('default');
rng(0);

imdb.imageDir = imgDir;

fid = fopen(fullfile(dbDir,isplit));
labels = [];
names = {};
classNames = {};
classId = 0;
newClass = 1;
while 1
    tline = fgetl(fid);
    if tline==-1
        break
    end   
    tline = strtok(tline,' ');
    
    if length(tline) == 7 && tline(1) >= '4' && tline(1) <= '6' 
        names = horzcat(names, [tline, '.avi']);
        labels = horzcat(labels, classId);
        newClass = 1;
        
    elseif newClass
        classId = classId + 1;
        newClass = 0;
        
        tline = fgetl(fid);
        if tline==-1
            break
        end
        tline = strtok(tline,' ');
        classNames = horzcat(classNames, tline);
        
    end
end
fclose(fid);

imdb.images.id = 1:numel(names);
imdb.images.name = names ;
imdb.images.set = zeros(1, numel(names)) ;
imdb.images.label = labels ;
imdb.classes.name = classNames;

% split training and testing sets
% trPerClass = 5;
for c=1:classId
    thisClass = find(labels==c);
    randClass = thisClass(randperm(length(thisClass)));
    trPerClass = fix((length(thisClass))/2);
    imdb.images.set(randClass(1:trPerClass)) = 1;
    imdb.images.set(randClass(trPerClass+1:end)) = 3;
    
end



end

