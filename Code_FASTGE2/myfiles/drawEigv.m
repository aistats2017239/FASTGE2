function drawEigv(x,nr,nc)
%% drawEigv plots the eigenvector x in a reshaped form
% Inputs:
%   x: eigenvector
%   nr: number of rows in the reshaped form
%   nc: number of columns in the reshaped form

%%
figure;
imagesc(reshape(x,nr,nc));
axis off;
set(gca,'position',[0 0 1 1],'units','normalized');
end