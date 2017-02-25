clear variables; close all;
addpath('./myfiles');
addpath('./ncutfiles');
addpath('./eigensolvers');

tt = tic;

%% Input the following three arguments:
% Image name
imgName = 'flower';
% Number of clusters of partitioning
numClusters = 2;
% Read in point from user input(0) or file(1)
ifFile = 1;
% Compute eigenproblems
% flag = 1: lobpcg for matrix pair (LG,LH);
%      = 2: lobpcg for matrix pair (K,M);
%      = 3: lobpcg for matrix pair (K,M) with preconditioning;
flag = 123;

%% Read in image and show it
[img,imgo,nr,nc] = readImage(imgName);
n = nr*nc;
% figure(2);image(imgo);axis off;
% set(gca,'position',[0 0 1 1],'units','normalized');
% pause;

%% Input the constraints points
numLabelSets = numClusters;
pts = readinPoints(ifFile,imgName,numLabelSets);

%% Compute W, d, LG and LH
% disp('start to compute W');
tw = tic;
[W,d] = computeWd(img);
ttw = toc(tw);
tgh = tic;
[LG,LH] = createConstrTL(pts,numClusters,nr,nc,W,d);
ttgh = toc(tgh);
clear W d;

%% Computing an eigenproblem
teig = tic;
tol = 1e-4;
maxit = 10000;
maxitIn = 4;
b = 0;
mu = 0.001;
[X1,lam1,resHist1,res1,relr1,...
 X2,lam2,resHist2,res2,relr2,...
 X3,lam3,resHist3,res3,relr3] = computeEigTL(LG,LH,numClusters,...
 tol,mu,maxit,maxitIn,b,flag);

tteig = toc(teig);

%% Plot the data
plotTL(X1,X2,X3,numClusters,img,imgo,pts,flag);
t_total = toc(tt);
ttt = [ttw ttgh tteig t_total-ttw-ttgh-tteig t_total];
% disp('     W        LGLH      eig       other     total');
% disp(ttt);
    