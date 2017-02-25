function drawObject(label,Img)
%% drawObject filters the object out of the background
% Inputs:
%   label: clustering results
%   Img: original image

%%
figure;
[nr,nc,~] = size(Img);
labelBack = label(1,1);
for i = 1:nr
    for j = 1:nc
        if(label(i,j)==labelBack)
            Img(i,j,1) = 0;
            Img(i,j,2) = 0;
            Img(i,j,3) = 0;
        end
    end
end
image(Img);
axis off;
set(gca,'position',[0 0 1 1],'units','normalized');
end