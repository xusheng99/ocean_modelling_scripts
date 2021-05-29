'''SPTDN算法'''
'''the Shortest Path of Time-Dependent Networks'''

import numpy as np

# 始末位置网格编号 #               # 网格编号从1开始计，python中数组和列表是从0开始计数的，所以涉及到数组（列表）中的定位要减一
start_day = -1                  # start_day = -1：不计初始时刻，找出所有start到end的路径；否则就是给定初始时刻。从0开始计数
start = 1                       # 起始网格编号
end = 6                         # 目标网格编号
N_nodes = 6                     # 网格数目

# 二维传输矩阵 #
inf = 10086                     # 给无限大定义一个值
# TEST1
mgraph1 = np.array([[0, 1, 12, inf, inf, inf],
          [inf, 0, 9, 3, inf, inf],
          [inf, inf, 0, inf, 5, inf],
          [inf, inf, 4, 0, 13, 15],
          [inf, inf, inf, inf, 0, 4],
          [inf, inf, inf, inf, inf, 0]])
mgraph2 = np.array([[0, 2, 10, inf, inf, inf],
           [inf, 0, inf, 6, 3, inf],
           [inf, 8, 0, inf, 4, inf],
           [inf, inf, inf, 0, 15, 13],
           [inf, inf, inf, inf, 0, 7],
           [inf, inf, inf, inf, inf, 0]])
mgraph3 = np.array([[0, inf, 7, inf, inf, inf],
           [4, 0, 11, 8, inf, inf],
           [inf, inf, 0, inf, 6, inf],
           [inf, inf, inf, 0, 11, 20],
           [inf, inf, inf, inf, 0, inf],
           [inf, inf, inf, inf, inf, 0]])
mgraph4 = np.array([[0, 4, 13, inf, inf, inf],
           [inf, 0, inf, 2, inf, inf],
           [inf, inf, 0, inf, 6, inf],
           [inf, inf, inf, 0, inf, 7],
           [inf, inf, inf, inf, 0, 8],
           [inf, inf, inf, inf, inf, 0]])
mgraph5 = np.array([[0,inf,inf,inf,inf,inf],
           [1,0,inf,7,inf,inf],
           [inf,4,0,inf,inf,inf],
           [inf,inf,8,0,10,inf],
           [inf,inf,inf,inf,0,4],
           [inf,inf,inf,inf,inf,0]])
mgraph = [mgraph1, mgraph2, mgraph3, mgraph4, mgraph5]      
'''
# TEST2
mgraph1 = np.array([[0, 1, inf, inf, inf, inf],
          [inf, 0, inf, inf, inf, inf],
          [inf, inf, 0, inf, inf, inf],
          [inf, inf, inf, 0, inf, inf],
          [inf, inf, inf, inf, 0, inf],
          [inf, inf, inf, inf, inf, 0]])
mgraph2 = np.array([[0, inf, inf, inf, inf, inf],
          [inf, 0, inf, 1, inf, inf],
          [inf, inf, 0, inf, inf, inf],
          [inf, inf, inf, 0, inf, inf],
          [inf, inf, inf, inf, 0, inf],
          [inf, inf, inf, inf, inf, 0]])
mgraph3 = np.array([[0, inf, inf, inf, inf, inf],
          [inf, 0, inf, inf, inf, inf],
          [inf, inf, 0, inf, inf, inf],
          [inf, inf, inf, 0, 1, inf],
          [inf, inf, inf, inf, 0, inf],
          [inf, inf, inf, inf, inf, 0]])
mgraph4 = np.array([[0, inf, inf, inf, inf, inf],
          [inf, 0, inf, inf, inf, inf],
          [inf, inf, 0, inf, inf, inf],
          [inf, inf, inf, 0, inf, inf],
          [inf, inf, 1, inf, 0, inf],
          [inf, inf, inf, inf, inf, 0]])
mgraph5 = np.array([[0, inf, inf, inf, inf, inf],
          [inf, 0, inf, inf, inf, inf],
          [inf, inf, 0, 1, inf, inf],
          [inf, inf, inf, 0, inf, inf],
          [inf, inf, inf, inf, 0, inf],
          [inf, inf, inf, inf, inf, 0]])
mgraph6 = np.array([[0, inf, inf, inf, inf, inf],
          [inf, 0, inf, inf, inf, inf],
          [inf, inf, 0, inf, inf, inf],
          [inf, inf, inf, 0, inf, 1],
          [inf, inf, inf, inf, 0, inf],
          [inf, inf, inf, inf, inf, 0]])
mgraph = [mgraph1, mgraph2, mgraph3, mgraph4, mgraph5, mgraph6]
'''
N_times = len(mgraph)           # 时间跨度

# 距离路径初始化 #
fit = np.ones((N_nodes, N_times))*inf                              # 从t时刻出发，i到终点end的距离，初始都为无限大
fit[end-1,:] = 0                                                    # 从t时刻出发，end到终点end的距离=0
pit = [[[end]for _ in range(N_times)]for _ in range(N_nodes)]      # 从t时刻出发，i到终点end的路径

# 开始寻找路径 #
vj = end-1                                                         # 后置节点vj的位置
lst_i = []
lst_j = []
for t in range(N_times):                                           # 在时间跨度中循环
    for vi in np.where((mgraph[t][:, vj] != inf) & (mgraph[t][:, vj] != 0))[0]:  # 找出t时刻vj的前置节点vi
        fit[vi, t] = mgraph[t][vi, vj]                             # 更新vi到end的距离
        pit[vi][t].insert(0, vi+1)                                 # 更新vi到end的路径
        # print(t, vi, fit[vi, t])
        lst_i.append(vi)                                           # 将vi添加到前置节点列表
        # 为什么以end为后置节点的要单独处理？
        # 因为：1.第一步时，所有路径都为无限大，所以不需要判断，可直接赋值；2.在更新路径时，作为最后的节点，因为缺少tt+1时刻的数据，虽然是end到end的距离，都为0，但时间对不上。
        for tt in range(t-1, -1, -1):                              # 逆向循环，从t时刻开始从后往前找
            lst_j = lst_i                                          # 上一步的前置节点变为这一步的后置节点
            lst_i = []                                             # 这一步的前置节点变为空集
            for vjj in lst_j:                                      # 在后置节点列表中循环，对于每一个后置节点vjj去寻找他的前置节点
                for vii in np.where((mgraph[tt][:, vjj] != inf) & (mgraph[tt][:, vjj] != 0))[0]:  # 找出tt时刻vjj的前置节点vii
                    if fit[vii, tt] > mgraph[tt][vii, vjj] + fit[vjj, tt+1]:                      # 判断，如果现在已知的tt时刻i到end的距离 > tt+1时刻j到end的距离 + tt时刻i到j距离，则更新；否则就停止停止运算
                        fit[vii, tt] = mgraph[tt][vii, vjj] + fit[vjj, tt+1]                      # 更新距离
                        pit[vii][tt] =[vii+1] + pit[vjj][tt+1]                                    # 更新路径
                        lst_i.append(vii)                                                         # 将找到的前置节点vii添加进列表

# 输出结果 #
# print(fit)
# print(pit)
SD = np.where(fit[start-1, :] < inf)[0]                            # 找出距离数组中，从i到end的距离不为无限大的位置
if start_day == -1:                                                # 输出全部路径
    for n in SD:
        # print(n)
        print('start day = {}: path ='.format(n), pit[start-1][n], 'dis =', fit[start-1, n], 'time =', len(pit[start-1][n]))
else:
    if start_day in SD:                                            # 输出指定初始时刻的路径
        print('start day = {}: path ='.format(start_day), pit[start - 1][start_day], 'dis =', fit[start - 1, start_day], 'time =', len(pit[start-1][start_day]))
    else:                                                          # 未找到路径
        print('start day = {}: No path for {} to {}!'.format(start_day, start, end))

'''
参考文献：
[1]谭国真,高文.时间依赖的网络中最小时间路径算法[J].计算机学报,2002(02):165-172.
'''