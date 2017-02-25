function [X1,lam1,resHist1,res1,relr1,...
          X2,lam2,resHist2,res2,relr2,...
          X3,lam3,resHist3,res3,relr3] = ...
          computeEigTL(LG,LH,k,tol,mu,maxit,maxitIn,b,flag)
%% computeEigTL solves LG x = lambda LH x in three different methods:
% (1) directly run LOBPCG
% (2) run LOBPCG after a spectral transformation
% (3) run LOBPCG with preconditioning after a spectral transformation
% Inputs:
%   LG: constrained graph Laplacian
%   LH: constrained graph Laplacian
%   k: LOBPCG block size
%   tol: tolerance of LOBPCG
%   mu: parameter of spectral transformation, used in method (2) and (3)
%   maxit: maximum number of LOBPCG iterations
%   maxitIn: maximum number of inner PCG iterations
%   b: extra LOBPCG block size
%   flag: indicating which method to be called
%         = 1, method (1)
%         = 2, method (2)
%         = 3, method (3)
%         It can be a combination, e.g. flag = 12 or 123
%         runing both method (1) and (2), or run all of them.
% Outputs: (for method (i), i = 1, 2 or 3)
%   Xi: eigenvectors corresponding to the k smallest eigenvalues
%   lami: k smallest lambda of LG x = lambda LH x
%   resHisti: residual history
%   resi: the absolute residual at the end of iterations
%   relri: the relative residual at the end of iterations

%% Initialization 
global gLG gLH gZ gMu gMaxitIn
gMu = mu;
n = size(LG,1);
gMaxitIn = maxitIn;
disp('start to compute the eigenvalue problem');
blSize = k + b ;
X0 = sprand(n,blSize,0.1);

%% (1) LOBPCG for (LG, LH)
if flag == 1 || flag == 12 || flag == 13 || flag == 123
    t1 = tic;
    [X1,lam1,~,~,resHist1] = lobpcg(X0,LG,LH,tol,maxit);
    res1 = resHist1(:,end);
    relr1 = zeros(blSize,1);
    for i = 1:blSize
        relr1(i) = res1(i)/(norm(LG*X1(:,i)) + lam1(i)*norm(LH*X1(:,i)));
    end
    disp(['LOBPCG(LG,LH) ' num2str(toc(t1)) ' seconds']);
else
    X1 = 0; lam1 = 0; resHist1 = 0; res1 = []; relr1 = [];
end

%% (2) LOBPCG for (K,M)
gLG=LG;
gLH=LH;
Z = ones(n,1);
gZ = Z/sqrt(n);

tolNull=[];
if ~isempty(tolNull)
    gZ = commonNullspace(gLG,gLH,tolNull);
end
if flag == 2 || flag == 12 || flag == 23 || flag == 123
    t2 = tic;
    [X2,sig2,~,~,resHist2] = lobpcg...
        (X0,-gLH,@regularMX,[],[],tol,maxit);
    lambdaOrigin = -1 ./ sig2 - mu;
    lam2 = lambdaOrigin;
    res2 = zeros(blSize,1);
    for i = 1:blSize
        res2(i)=norm(LG*X2(:,i) - LH*X2(:,i)*lam2(i));
    end
    relr2 = zeros(blSize,1);
    for i = 1:blSize
        relr2(i) = res2(i)/(norm(LG*X2(:,i)) + lam2(i)*norm(LH*X2(:,i)));
    end
    disp(['LOBPCG(K,M) ' num2str(toc(t2)) ' seconds']);
else
    X2 = 0; lam2 = 0; resHist2 = 0; res2 = []; relr2 = [];
end

%% (3) LOBPCG for (K,M) with preconditioner
if flag == 3 || flag == 13 || flag == 23 || flag == 123
    t3 = tic;
    if ~isempty(tolNull)
        gZ = commonNullspace(gLG,gLH,tolNull);
    end
    [X3,sig3,~,~,resHist3] = lobpcg...
        (X0,-gLH,@regularMX,@myLinCG,[],tol,maxit);
    lam3 = -1 ./ sig3 - mu;
    res3 = zeros(blSize,1);
    for i = 1:blSize
        res3(i) = norm(LG*X3(:,i) - LH*X3(:,i)*lam3(i));
    end
    relr3 = zeros(blSize,1);
    for i = 1:blSize
        relr3(i) = res3(i)/(norm(LG*X3(:,i)) + lam3(i)*norm(LH*X3(:,i)));
    end
    disp(['LOBPCG(K,M) with precond. '...
        num2str(toc(t3)) ' seconds']);
else
    X3 = 0; lam3 = 0; resHist3 = 0; res3 = []; relr3 = [];
end
end
