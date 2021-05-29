ax = gca;
mycmap = colormap(ax); 
save('MyColormaps','mycmap');
%保存colorar
ax = gca;
load('MyColormaps','mycmap')
colormap(ax,mycmap)
%使用colorbar