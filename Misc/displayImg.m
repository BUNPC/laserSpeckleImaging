function displayImg(img,perc,imgTitle,cmap)
imagesc(img);
colormap(cmap)
xticklabels([]);
yticklabels([]);
colorbar
caxis([prctile(img(:),perc(1)),prctile(img(:),perc(2))]);
title(imgTitle);
axis image;
end