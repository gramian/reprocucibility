
cmap =colormap;
lcmap = size(cmap,1);

 val  = objvalD(sub2ind(size(objvalD),I,J));
valmax = max(val);
valmap = round(val./(valmax) * (lcmap-1))+1;

figure
hold on
for i = 1:length(x)
    plot(x(i),y(i),'.','markeredgecolor',cmap(valmap(i),:))
end

    
