import numpy as np
import cv2
#读取一张图片
size = (432,288)
print(size)
#完成写入对象的创建，第一个参数是合成之后的视频的名称，第二个参数是可以使用的编码器，第三个参数是帧率即每秒钟展示多少张图片，第四个参数是图片大小信息
videowrite = cv2.VideoWriter(r'F:\test.mp4',-1,20,size)#20是帧数，size是图片尺寸
img_array=[]
for filename in [r'F:\Picture\{0}.png'.format(i) for i in range(600)]:
 img = cv2.imread(filename)
 if img is None:
  print(filename + " is error!")
  continue
 img_array.append(img)
for i in range(600):
 videowrite.write(img_array[i])
print('end!')