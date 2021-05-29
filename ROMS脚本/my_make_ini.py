# 困难: 
# 1. z2sigma的函数写不出来, 或者调用matlab中的zlevs?
# 2. 数据源问题, 很难下载到高时间分辨率和空间分辨率的数据, HYCOM极难下载, SODA又缺失16年之后的数据.
# 3. 对于grid.nc中没有值的位置, 如何填充一个还近似合理的值?
# 4. ini中不仅有温盐流等数据, 还有一些其他的小变量, 这些变量如何设置?

import os
from netCDF4 import Dataset
import numpy as np
import math
import shutil


