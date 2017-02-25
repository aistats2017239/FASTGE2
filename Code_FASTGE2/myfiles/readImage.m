function [img,imgo,nr,nc] = readImage(imgName)
%% readImage reads the image into the matrix
% Inputs:
%   imgName: a string of image name
% Outputs:
%   img: a grey image converted from the original image
%   imgo: original image in matrix form
%   nr: number of rows of the image
%   nc: number of columns of the image

%%
img = imread(strcat('./data/',imgName,'.jpg'));
imgo = img;
[nr,nc,nb] = size(img);
if (nb > 1)
    img = double(rgb2gray(img));
else
    img = double(img);
end
figure(1);clf;imagesc(img);colormap(gray);
axis off;set(gca,'position',[0 0 1 1],'units','normalized');
end