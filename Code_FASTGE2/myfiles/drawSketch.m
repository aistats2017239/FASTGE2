function drawSketch(label,Img,pts,k)
%% drawSketch draws the boundary between different clusters
% Inputs:
%   label: clustering results
%   Img: original image
%   pts: constraint points
%   k: number of clusters

%%
colorSet = {'b','g','y','c','m'};
figure;
bw = edge(label,0.01);
J=showmask(Img,imdilate(bw,ones(4,4)));
imagesc(J);colormap(gray);
axis off;
set(gca,'position',[0 0 1 1],'units','normalized');
hold on;
for i = 1:k
    scatter(pts{i}(:,1),pts{i}(:,2),100,colorSet{i},'filled');
end
hold off;
end