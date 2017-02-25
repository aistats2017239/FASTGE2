function [LG,LH] = createConstrTL(pts,numClusters,nr,nc,W,d)
%% createConstrTL incorporates the constraints into two Laplacian LG and LH
% Inputs:
%   pts: constraint points
%   numClusters: number of clusters to be computed
%   nr: number of rows in the image
%   nc: number of columns in the image
%   W: affinity matrix
%   d: diagonal of the degree matrix
% Outputs: (for method (i), i = 1, 2 or 3)
%   LG: Laplacian with Must-Link constraints
%   LH: Laplacian with Cannot-Link constraints

%% Initialization
n = nr*nc;
m = zeros(numClusters,1);
v = cell(numClusters);
dmin = min(d);
dmax = max(d);
dmm = dmin * dmax;

%% Construct must-link constraint WM
Wm = sparse(n,n);
for k = 1:numClusters
    m(k) = size(pts{k},1);
    v{k} = zeros(m(k),1);
    p = round(pts{k});
    for i = 1:m(k)
        v{k}(i) = (p(i,1)-1)*nr+p(i,2);
        for j = 1:i-1
            Wm(v{k}(j),v{k}(i)) = d(v{k}(j))*d(v{k}(i))/dmm;
        end
    end
end
WM = Wm + Wm';

%% Construct cannot-link constraint WC
Wc = sparse(n,n);
for k = 1:numClusters
    for kk = k+1:numClusters
        for i = 1:m(k)
            for j = 1:m(kk)
                Wc(v{k}(i),v{kk}(j)) = d(v{k}(i))*d(v{kk}(j))/dmm;
            end
        end
    end
end
WC = Wc + Wc';
dc = sum(abs(WC),2);
volc = sum(dc);

%% Form the Laplacians
WG = W + 3*WM;
dG = sum(WG,2);
DG = spdiags(dG,0,n,n);
LG = DG - WG;

WH = 3*WC + dc*dc'/(volc*n);
dH = sum(WH,2);
DH = spdiags(dH,0,n,n);
LH = DH-WH;
end