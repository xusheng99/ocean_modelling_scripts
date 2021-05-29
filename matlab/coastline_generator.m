%lon1,lon2,lat1,lat2分别指目标区域的左右下上边界的经纬度，以小数或整数格式输入
%quality参数可选'f','h','i','l','c'
%最后一个参数为生成岸线文件的名字

lonmin =  105;   % Minimum longitude [degree east]
lonmax =  140;   % Maximum longitude [degree east]
latmin =   5;   % Minimum latitudeF  [degree north]
latmax =   41;   % Maximum latitude  [degree north]

coastline_gen(lonmin,lonmax,latmin,latmax,'f','all_china_sea_f.mat')



function [] = coastline_gen(lon1,lon2,lat1,lat2,quality,name)

addpath 'C:\Users\XUSHENG\software\m_map'    %将m_map的文件夹添加到路径
m_proj('Equidistant Cylindrical','lat',[lat1 lat2],'long',[lon1 lon2]);
if quality == 'c'
    m_gshhs_c('save','topodata_c') %将所需岸线保存为topodata.mat文件
end
if quality == 'l'
    m_gshhs_l('save','topodata_l')
end
if quality == 'i'
    m_gshhs_i('save','topodata_i')
end
if quality == 'h'
    m_gshhs_h('save','topodata_h')
end
if quality == 'f'
    m_gshhs_f('save','topodata_f')
end

oldname = ['topodata_',quality,'.mat']
movefile(oldname,name)

load(name);  %连续岸线
plot(ncst(:,1),ncst(:,2));    %两个同长度一维数组绘制连线图 scatter为散点，plot则将这些点顺序连接
axis equal
axis([lonmin,lonmax,latmin,latmax]);    %岸线数据范围是全球，我们只画目标区域的

end

