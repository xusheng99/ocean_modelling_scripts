## DIJKSTRA

**问题是什么? **
给定图, 指定两个节点, 求路径P, 使得P经过的边的权和最小

**HINTS?**
如果A-B-C-D-E-F-G这样的一条路径已经被确定, 那么其中的任意一条子路径, 例如A-B-C, 或 B-C-D-E-F也是A到C/B到F的最短路径.
核心定理如上.

问题转化成:
找到一条路径, 所有的点都会被不遗漏的经过. 而且边权和最小.(即**邮递员问题**)
这样的转化似乎是把原问题复杂化了,  加上了新的限制条件(需要经过所有点)
但实际上在添加了新的限制条件的同时也削弱了原问题的限制, 即从A出发, 但并不限制路径在何处结束, 而原问题是要求A出发至G结束的.

只要找到了邮递员的最佳路径, 那么可以容易的从中取出任意两点之间的最短路径.
那么, 新问题成功被**转化**成新问题: 邮递员问题

**邮递员问题作何解答呢?**

邮递员问题相比与原问题要更容易的地方在于: 邮递员并没有明确的目标节点, 只要遍历所有节点就可以了.
所以, 我们可以像蜗牛用触手探路一般, 在已确定最短路径的节点的集合周围, 触摸找出所有的临近节点.
再从这些临近节点中找到节点集合最近的那一个, 把它作为新的成员加进来, 吃掉它. 
可以肯定, 从初始节点到这个最新的节点, 集合中确定的路径**仍然**是最短的.

如此重复, 直到所有的点都被我们这只小蜗牛的触手摸过并吃掉, 那么邮递员的任务就被我们完成了.

完成邮递员任务后, 原问题也作为简单的推论自然得解.



**最快(短)路径和最可能路径的关系**

首先, Dijkstra算法计算的是边权和最小的路径. 若边权代表节点之间的通行时间,  则D算法求出的是最快的从指定节点出发遍历所有节点的路径(容易从中提取出指定两点之间的路径). 

下面约定一个记号: e(A,B,G) = value  意为在图G中, A节点和B节点之间存在由A指向B的边, 且边权为value.

在我们构建的流体网络中, 边权代表的是两个节点之间输运粒子的数目. e(A,B,G1)=9, 意为在A节点释放的若干粒子中, 有9个被输运到了B节点. 将每个节点的所有边的度都除以单个节点的释放粒子总数, 得到的新的邻接矩阵:

G2 = G1/(release\_particle\_counts)

G2和G1的底图是同构的, 唯一改变的只有边权.

新的邻接矩阵的度可做如下解释, e(A,B,G2) = 9/300 =0.03, 代表从A释放的所有粒子中, 有%3被输运到了B节点. 或者另一种理解为, A节点释放的粒子有0.03的**概率**会到B节点.

对于G2矩阵, 使用dijkstra算法计算出来的是**概率和**最小的路径. 概率相加并没有什么意义, 况且还是概率和**最小**的那条路径. 这和我们的目标: 最可能路径 是背道而驰的. 因此, 要么采用构建一个新的矩阵G3, G3和G2同构, 且G3下的**边权和最小**路径恰为G2下的**边权积最大**路径, 要么就采用新的算法, 使之能适应我们的要求.

不妨试试构建一个新矩阵的方法先.

使用负自然对数可以满足这样的要求. dijkstra算法不支持负边权图的计算. 底数为e大于1, 保证对于所有的边权取对数后得到的都是负值, 再取相反数, 全为正值.






