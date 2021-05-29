file='uv.dat';
load(file);

lon=uv(:,1);
lat=uv(:,2);
u=uv(:,3);
v=uv(:,4);

for i = 1:20:length(lon);
        quiver(lat(i,1),lon(i,1),u(i,1),v(i,1));
        hold on;
end;

