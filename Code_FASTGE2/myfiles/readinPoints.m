function pts = readinPoints(ifFile,imgName,numLabelSets)
%% read the points from user's input or a file
% Inputs:
%   ifFile: specifying the way to input the points
%           = 0, user inputs the points
%           = 1, read points from a file
%   imgName: a string of image name
%   numLabelSets: number of sets/groups of points
% Outputs:
%   pts: constraint points

%%
colorSet = {'b','g','y','c','m'};
if (ifFile == 1)
    % Read 'pts' in constraint points
    load(strcat('./data/pts_', imgName, '.mat'));
    figure(1);
    hold on;
    for i = 1:numLabelSets
        scatter(pts{i}(:,1),pts{i}(:,2),100,colorSet{i},'filled');
    end
    hold off;
else
    figure(1);
    hold on;
    for i = 1:numLabelSets
        [x, y] = ginput;
        pts{i} = [x,y];
        scatter(pts{i}(:,1),pts{i}(:,2),100,colorSet{i},'filled');
    end
    hold off;
end
end