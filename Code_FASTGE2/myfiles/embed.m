function [c,Y] = embed(X,k)
%% embed renormalizes X into Y and clusters the rows of Y
% Inputs:
%   X: the computed eigenvector
%   k: number of clusters
% Outputs:
%   c: labels of clusters
%   Y: renormalized eigenvector

%%
n = size(X,1);
if k <= 2
    c = X(:,1) > 0;
    c = c+1;
    XX = X(:,1:k);
    for i = 1:k
        XX(:,k) = XX(:,k)/norm(XX(:,k));
    end
    Y = zeros(n,k);
    for i = 1:n
        Y(i,:) = XX(i,:)/norm(XX(i,:));
    end
else
    k = k-1;
    XX = X(:,1:k);
    for i = 1:k
        XX(:,k) = XX(:,k)/norm(XX(:,k));
    end
    Y = zeros(n,k);
    for i = 1:n
        Y(i,:) = XX(i,:)/norm(XX(i,:));
    end
    c = kmeans(Y,k+1);
end
end
