function [W,d] = computeWd(Img)
%% computeWd computes the affinity matrix and the degree matrix of Img
% Inputs:
%   Img: original image data
% Outputs:
%   W: affinity matrix
%   d: degree matrix

%% Compute semilarity matrix
W = ICgraph(Img);
W = sparsifyc(W,1e-6);
% check for matrix symmetry
% if max(max(abs(W-W'))) > 1e-10 %voir (-12) 
%     %disp(max(max(abs(W-W'))));
%     error('W not symmetric');
% end

%% Compute the diagonal of degree matrix
d = sum(abs(W),2);
end
