function MX=regularMX(X)
%% regularMX forms the matrix-vector multiplication M * X
% where M = LG + mu*LH + Z*Z'
%  Inputs:
%    X: the projected matrix in Rayleigh Ritz method
%  Outputs:
%    MX: the matrix-vector multiplication

%%
global gLG gLH gZ gMu
MX = gLG * X + gMu * (gLH * X) + gZ * (gZ' * X);
end