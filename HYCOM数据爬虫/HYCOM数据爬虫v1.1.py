import numpy as np
import os
import urllib
import urllib3
import requests
import arrow
import time
import random
import sys

'''
v1.1 更新日志
加入自动检测缺失文件并重新下载的功能
加入命令行中指定年份和存储文件夹的功能, 避免反复修改源代码
'''

targetyear = 2019 #目标年份
workdir = '' #工作目录 
# sys.argv[] 而不是 sys.argv()!!!

#给定开始和结束的序号 加入跳过已下载功能后这里就不再需要修改了
start_num = 0  #起始 闭区间
end_num = 365  #结束 开区间
    
os.chdir(workdir)
pwd = os.getcwd()
print('now working at this dir: ',workdir,'please check, all files will be stored to this folder')

def alldatelist(year):      
    year = str(year)
    def isLeapYear(years):
        '''
        通过判断闰年，获取年份years下一年的总天数
        :param years: 年份，int
        :return:days_sum，一年的总天数
        '''
        # 断言：年份不为整数时，抛出异常。
        assert isinstance(years, int), "请输入整数年，如 2018"
        
        if ((years % 4 == 0 and years % 100 != 0) or (years % 400 == 0)):  # 判断是否是闰年
            days_sum = 366
            return days_sum
        else:
            days_sum = 365
            return days_sum
    def getAllDayPerYear(years):
        '''
        获取一年的所有日期
        :param years:年份
        :return:全部日期列表N
        '''
        start_date = '%s-1-1' % years
        a = 0
        all_date_list = []
        days_sum = isLeapYear(int(years))
        print()
        while a < days_sum:
            b = arrow.get(start_date).shift(days=a).format("YYYY-MM-DD")
            a += 1
            all_date_list.append(b)
        return all_date_list
    if __name__ == '__main__':
        all_date_list = getAllDayPerYear(year)
        return all_date_list

def date(x):    #某个year的第x位日子，x+1，天数+5
    #以2019-01-01为第一个日子，x每加1，日期加五天，要求以"year-month-day"格式输出
    #返回值为形如"2020-01-26"的格式
    interval = 1
    date = alldatelist(targetyear)[interval*x]   
    return date

def ncname(x):
    #每个文件名字不能相同，否则后面下载的会覆盖前面的
    name = str(4*x+1000)+'_'+date(x)+'.nc4'
    return name

def urlgen(x):  #按hycom格式生成所需要的下载链接
    # 参数直接在链接里改
    prefix = 'https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0?var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v&north=26&west=105&east=125&south=5&disableProjSubset=on&horizStride=1&time='
    datepart = date(x)
    suffix = 'T00%3A00%3A00Z&vertStride=1&addLatLon=true&accept=netcdf4'
    return(prefix+datepart+suffix)

brokenfile_list = []
def download(i):  #下载目标日期序列中第i个日期的nc文件
    print('download the',i+1,'th file')
    print('data:',date(i))

    # 向HYCOM网站请求文件
    url = urlgen(i)
    r = requests.get(url)
    print('requests done')

    #写入文件
    with open(ncname(i),'wb') as f:
        f.write(r.content)
        f.close()

    # 判断文件是否损坏
    if os.path.getsize(ncname(i)) < 10000:
        os.remove(ncname(i))
        brokenfile_list.append(i) #标注损坏文件, 稍后重新下载
        print('broken file: size too small!')
        print('the broken file has been removed!')
    else:
        print('download completed, now wait for about 2 minutes to cheat the server')

    # 暂停5分钟, 防止封IP  经过测试, 这个暂停时间很稳, 基本没断连过.
    time.sleep(300+((random.randint(1,120)-60)/2))
    print('stopped for about 5 minutes,now continue')
    print('\n')

if __name__ == "__main__":
    try:
        filelist = os.listdir(workdir)
        for i in range(start_num,end_num):
            if ncname(i) in filelist:
                print(ncname(i),': this file already been downloaded! It will be skipped!')
            else: 
                download(i)
        print(brokenfile_list)

        if len(brokenfile_list) != 0:
            for j in brokenfile_list:
                download(j)
        else:
            print('mission this time is fullly completed with no errors!')
            
        #报个数
        fl = os.listdir()
        fllen = len(fl)
        print('downloaded',str(fllen),'files')
    except ConnectionResetError or WindowsError or Sizeerror:
        print('failed')
    finally:
        print('fucked hycom.org')
