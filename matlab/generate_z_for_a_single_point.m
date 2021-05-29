clc;
clear;

h = ncread('his_00161.nc','h')
lonlat = load('lonlat.txt')
lonlat = lonlat+1

for i = 1:36
    lon = lonlat(i,1)
    lat = lonlat(i,2)
    
   
    
    if lon == 190
        lonreal = 115.5
    end
    if lon == 192
        lonreal = 115.6
    end
    if lon == 194
        lonreal = 115.7
    end
    if lon == 195
        lonreal = 115.8
    end
    if lon == 197
        lonreal = 115.9
    end
    if lon == 199
        lonreal = 116.0
    end
    if lat == 280
        latreal = 20.5
    end
    if lat == 282
        latreal = 20.6
    end
    if lat == 284
        latreal = 20.7
    end
    if lat == 285
        latreal = 20.8
    end
    if lat == 287
        latreal = 20.9
    end
    if lat == 289
        latreal = 21.0
    end

zname = [num2str(lonreal,'%6.1f'),'E_',num2str(latreal','%6.1f'),'N.txt'] 
hp = h(lon,lat);
z = zlevs(hp,0,4.5,1.5,5,40,'r',2);
z = -z;
zz = 1:40;
for k = 1:40
    zz(k) = z(41-k);
end
z = zz

outfile = fopen(zname,'w');
for line = 1:40
    fprintf(outfile,'%15.10f',z(line));
    fprintf(outfile,'\n');
end
fclose(outfile);


end

    

   