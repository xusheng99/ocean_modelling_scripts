import numpy as np
import os
import urllib
import urllib3
import requests
import arrow
import time
import random

#给定开始和结束的序号
start_num = 0  #起始 闭区间
end_num = 254  #结束 开区间

os.chdir('C:\\Users\\XUSHENG\\Desktop\\hycom_DL\\')



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
            # print(years, "是闰年")
            days_sum = 366
            return days_sum
        else:
            # print(years, '不是闰年')
            days_sum = 365
            return days_sum
    def getAllDayPerYear(years):
        '''
        获取一年的所有日期
        :param years:年份
        :return:全部日期列表
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
        # print(all_date_list)
        return all_date_list
    if __name__ == '__main__':
        # years = "2001"
        # years = int(years)
        # # 通过判断闰年，获取一年的总天数
        # days_sum = isLeapYear(years)
        # 获取一年的所有日期
        all_date_list = getAllDayPerYear(year)
        return all_date_list

def date(x):    #某个year的第x位日子，x+1，天数+5
    #以2019-01-01为第一个日子，x每加1，日期加五天，要求以"year-month-day"格式输出
    #返回值为形如"2020-01-26"的格式
    interval = 1
    targetyear = 2020  #####  !!!!!  #####  !!!!!  年份在这里写
    date = alldatelist(targetyear)[interval*x]   
    return date

def ncname(x):
    #每个文件名字不能相同，否则后面下载的会覆盖前面的
    name = str(x+1000)+'_'+date(x)+'.nc4'
    return name

def urlgen(x):  #按hycom格式生成所需要的下载链接
    # 参数直接在链接里改
    prefix = 'https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0?var=surf_el&var=salinity&var=water_temp&var=water_u&var=water_v&north=41&west=117&east=132&south=24&disableProjSubset=on&horizStride=3&time='
    datepart = date(x)
    suffix = 'T00%3A00%3A00Z&vertStride=2&addLatLon=true&accept=netcdf4'
    return(prefix+datepart+suffix)

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
        print('broken file: size too small!')
        print('the broken file has been removed!')
    else:
        print('download completed, now wait for about 2 minutes to cheat the server')

    # 暂停5分钟, 防止封IP
    time.sleep(300+((random.randint(1,120)-60)/2))
    print('stopped for about 2 minutes,now continue')
    print('\n')

if __name__ == "__main__":
    try:
        for i in range(start_num,end_num):
            download(i)
        pass
        #报个数
        fl = os.listdir()
        fllen = len(fl)
        print('downloaded',str(fllen),'files')
    except ConnectionResetError or WindowsError or Sizeerror:
        print('failed')
    finally:
        print('fucked hycom.org')

