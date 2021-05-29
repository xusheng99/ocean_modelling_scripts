load('uv.dat');
lat=uv(:,1);
lon=uv(:,2);
u=uv(:,3);
v=uv(:,4);
speed=uv(:,5);

outfile=fopen('speed.txt','w');
for i = 1:length(lon);
    fprintf(outfile,'%15.10f %15.10f %20.10f',lon(i,1),lat(i,1),speed(i,1));
    fprintf(outfile,'\n');
end;
fclose(outfile);

outfile=fopen('u_v.txt','w');
for i = 1:length(lon);
    fprintf(outfile,'%15.10f %15.10f %20.10f %20.10f',lon(i,1),lat(i,1),u(i,1),v(i,1));
    fprintf(outfile,'\n');
end;
fclose(outfile);

