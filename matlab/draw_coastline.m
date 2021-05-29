function draw_coastline()

coast = load('./llbounds.bln');
index = find(coast(:,2) == 1);

for i = 1:(length(index)-1)
    polygon = coast((index(i)+1:index(i+1)-1),:);
    plot(polygon(:,1),polygon(:,2));
    hold on
end


