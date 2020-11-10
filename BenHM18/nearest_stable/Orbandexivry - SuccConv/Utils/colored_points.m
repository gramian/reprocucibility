cmap =colormap(jet(256));
% cmap(1,:) = [1 1 1];
lcmap = size(cmap,1);

val  = objvalD(sub2ind(size(objvalD),I,J,K));
valmax = max(val);
valmap = round((val-minnorm)./(valmax-minnorm) * (lcmap-1))+1;

figure
hold on
for i = 1:length(x)
    plot3(x(i),y(i),z(i),'o','MarkerSize',6,'markerfacecolor',cmap(valmap(i),:) ,'markeredgecolor',cmap(valmap(i),:))
end

    
