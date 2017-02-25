clear variables; close all;
addpath('./myfiles');
%addpath('./ncutfiles');
addpath('./eigensolvers');

%% Initialization
global gLG gLH gZ gMu
numEigen = 2;
mu = 0.001;

%% Create data points
ns = 100;
offset = 10;

S1 = rand(ns,2);
S1(:,1) = S1(:,1)*5;
S1(:,2) = S1(:,2)*10;

S2 = rand(ns,2);
S2(:,1) = S2(:,1)*5+offset;
S2(:,2) = S2(:,2)*10;

S3 = rand(ns,2);
S3(:,1) = S3(:,1)*5+offset*2;
S3(:,2) = S3(:,2)*10;

S = [S1;S2;S3];
n = size(S,1);

%% Compute affinity matrix and Laplacian
sigma = 0.5;
W = zeros(n,n);
for i = 1:n
    for j = 1:n
        W(i,j) = exp(-sqrt(sum((S(i,:)-S(j,:)).^2))/(2*sigma^2));
    end
end
d = sum(W,2);
D = diag(d);
L = diag(d)-W;

%% Create the constraint matrix
ns1 = 1;
ns2 = 2;
ns3 = 1;
ptInd1 = ceil(rand(ns1,1)*ns);
ptInd2 = ceil(rand(ns2,1)*ns)+ns;
ptInd3 = ceil(rand(ns3,1)*ns)+2*ns;

pt1 = [ptInd1;ptInd3];
pt2 = ptInd2;
m1 = size(pt1,1);
m2 = size(pt2,1);
dmin = min(d);
dmax = max(d);
dmm = dmin*dmax;

% Must-link constraints
Wm1 = zeros(n,n);
for i = 1:m1
    for j = i+1:m1
        Wm1(pt1(i),pt1(j)) = d(pt1(i))*d(pt1(j))/dmm;
        Wm1(pt1(j),pt1(i)) = d(pt1(i))*d(pt1(j))/dmm;
    end
end

Wm2 = zeros(n,n);
for i = 1:m2
    for j = i+1:m2
        Wm2(pt2(i),pt2(j)) = d(pt2(i))*d(pt2(j))/dmm;
        Wm2(pt2(j),pt2(i)) = d(pt2(i))*d(pt2(j))/dmm;
    end
end

% Cannot-link constraints
Wc = zeros(n,n);
for i = 1:m1
    for j = 1:m2
        Wc(pt1(i),pt2(j)) = d(pt1(i))*d(pt2(j))/dmm;
        Wc(pt2(j),pt1(i)) = d(pt1(i))*d(pt2(j))/dmm;
    end
end

%% Form LG and LH
alpha = 1;
beta = 1;
WG = W + alpha*(Wm1+Wm2);
dG = sum(WG,2);
LG = diag(dG) - WG;

dc = sum(Wc,2);
mc = sum(dc);
WH = beta*Wc + dc*dc'/(2*mc*n);
dH = sum(WH,2);
LH = diag(dH) - WH;

%% Solve the eigenvalue problem by direct calling LOBPCG on (LG,LH)
maxit = 10000;
tol = 1e-4;
X0 = rand(n,numEigen);
[X1,lam1,~,~,resHist1] = lobpcg(X0,LG,LH,tol,maxit);

%% Solve the eigenvalue problem by applying LOBPCG on (K,M)
% where K = -LH and M = LG + mu*LH + Z*Z'
gLG = LG;
gLH = LH;
Z = ones(n,1);
gZ = Z/sqrt(n);
gMu = mu;

[X2,lam2,~,~,resHist2] = lobpcg...
    (X0,-gLH,@regularMX,[],[],tol,maxit);
lamOrg = -1./lam2 - mu;
lam2 = lamOrg;
res2 = zeros(numEigen,1);
for i = 1:numEigen
    res2(i)=norm(LG*X2(:,i)-LH*X2(:,i)*lam2(i));
end

%% plot data
area = 50;
lsize = 30;

% original data with constraints
figure;
scatter(S(:,1),S(:,2),area,'filled');
hold on;
scatter(S(pt1,1),S(pt1,2),area*3,'r','LineWidth',3);
scatter(S(pt2,1),S(pt2,2),area*3,'g','LineWidth',3);
hold off;
xlabel('x','FontSize',lsize);
ylabel('y','FontSize',lsize);
set(gca,'fontsize',lsize);
xlim([0,25]);

% x1 returned by lobpcg(L_G,L_H)
figure;
scatter((1:n)',X1(:,1),area,'filled');
xlabel('i','FontSize',lsize);
ylabel('x_1(i)','FontSize',lsize);
set(gca,'fontsize',lsize);
xlim([0,300]);

% x1 returned by lobpcg(K,M)
figure;
scatter((1:n)',X2(:,1),area,'filled');
xlabel('i','FontSize',lsize);
ylabel('x_1(i)','FontSize',lsize);
set(gca,'fontsize',lsize);
xlim([0,300]);

% clustering by x1 of matrix pencil (K,M)
figure;
scatter(S1(:,1),S1(:,2),area,'filled','r');
hold on;
scatter(S2(:,1),S2(:,2),area,'filled','g');
scatter(S3(:,1),S3(:,2),area,'filled','r');
hold off;
xlabel('x','FontSize',lsize);
ylabel('y','FontSize',lsize);
set(gca,'fontsize',lsize);
xlim([0,25]);