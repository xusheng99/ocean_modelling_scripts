from netCDF4 import Dataset
import numpy as np 
import os 
from matplotlib import pyplot as plt

'''
netcdf roms_bulk {
dimensions:
	longitude = 53 ;
	latitude = 47 ;
	time = 7201 ;
variables:
	float time(time) ;
		time:long_name = "Time" ;
		time:units = "seconds since 2000-01-01 00:00:00.0" ;
	float longitude(longitude) ;
		longitude:long_name = "longitude" ;
		longitude:units = "degrees_east" ;
	float latitude(latitude) ;
		latitude:long_name = "latitude" ;
		latitude:units = "degrees_north" ;
	float Uwind(time, latitude, longitude) ;
		Uwind:long_name = "surface u-wind component" ;
		Uwind:units = "meter second-1" ;
		Uwind:time = "time" ;
		Uwind:coordinates = "longitude latitude" ;
	float Vwind(time, latitude, longitude) ;
		Vwind:long_name = "surface v-wind component" ;
		Vwind:units = "meter second-1" ;
		Vwind:time = "time" ;
		Vwind:coordinates = "longitude latitude" ;
	float Tair(time, latitude, longitude) ;
		Tair:long_name = "surface air temperature" ;
		Tair:units = "Celsius" ;
		Tair:time = "time" ;
		Tair:coordinates = "longitude latitude" ;
	float Pair(time, latitude, longitude) ;
		Pair:long_name = "surface air pressure" ;
		Pair:units = "Pa" ;
		Pair:time = "time" ;
		Pair:coordinates = "longitude latitude" ;
	float Qair(time, latitude, longitude) ;
		Qair:long_name = "surface air relative humidity" ;
		Qair:units = "percentage" ;
		Qair:time = "time" ;
		Qair:coordinates = "longitude latitude" ;
	float cloud(time, latitude, longitude) ;
		cloud:long_name = "cloud fraction" ;
		cloud:units = "" ;
		cloud:time = "time" ;
		cloud:coordinates = "longitude latitude" ;
	float rain(time, latitude, longitude) ;
		rain:long_name = "rain fall rate" ;
		rain:units = "kilogram meter-2 second-1" ;
		rain:time = "time" ;
		rain:coordinates = "longitude latitude" ;
	float swrad(time, latitude, longitude) ;
		swrad:long_name = "solar shortwave radiation flux" ;
		swrad:units = "watt meter-2" ;
		swrad:time = "time" ;
		swrad:coordinates = "longitude latitude" ;
	float lwrad_down(time, latitude, longitude) ;
		lwrad_down:long_name = "downwelling longwave radiation flux" ;
		lwrad_down:units = "watt meter-2" ;
		lwrad_down:time = "time" ;
		lwrad_down:coordinates = "longitude latitude" ;
}
'''

def ncread(ncfile,varname):
    dataset = Dataset(ncfile)
    a = dataset.variables[varname][:]
    return a

def ncwrite(matrix,ncfile,varname):
    content = Dataset(ncfile,'r+') #务必r+
    content.variables[varname][:] = matrix
    print(ncfile+'\'s variable:',varname,'\n','has been written!')
    content.close()

