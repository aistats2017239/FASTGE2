function plotTL(X1,X2,X3,k,img,imgo,pts,flag)
%% plotTL plots the image segmentation based on the computed eigenvectors
% from different methods:
% (1) directly run LOBPCG
% (2) run LOBPCG after a spectral transformation
% (3) run LOBPCG with preconditioning after a spectral transformation
% Inputs:
%   X1: the computed eigenvector by method (1)
%   X2: the computed eigenvector by method (2)
%   X3: the computed eigenvector by method (3)
%   k: number of clusters
%   img: the grey image converted from the original image
%   imgo: the original image
%   pts: contraint points
%   flag: indicating which method to be called
%         = 1, method (1)
%         = 2, method (2)
%         = 3, method (3)
%         It can be a combination, e.g. flag = 12 or 123
%         runing both method (1) and (2), or run all of them.

%%
[nr,nc,~] = size(img);
if flag == 1 || flag == 12 || flag == 13 || flag == 123
    [c1,Y1] = embed(X1,k);
    label1 = reshape(c1,nr,nc);
    drawSketch(label1,img,pts,k);
    drawObject(label1,imgo);
    drawEigv(X1(:,1),nr,nc);
    drawEigv(Y1(:,1),nr,nc);
end

if flag == 2 || flag == 12 || flag == 23 || flag == 123
    [c2,Y2] = embed(X2,k);
    label2 = reshape(c2,nr,nc);
    drawSketch(label2,img,pts,k);
    drawObject(label2,imgo);
    drawEigv(X2(:,1),nr,nc);
    drawEigv(Y2(:,1),nr,nc);
end

if flag == 3 || flag == 13 || flag == 23 || flag == 123
    [c3,Y3] = embed(X3,k);
    label3 = reshape(c3,nr,nc);
    drawSketch(label3,img,pts,k);
    drawObject(label3,imgo);
    drawEigv(X3(:,1),nr,nc);
    drawEigv(Y3(:,1),nr,nc);
end
end
