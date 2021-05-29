clear;clc;

addpath c:/Users/XUSHENG/Desktop/scripts/flow_networks/
%% 载入数据
A = load('Adjacency_Matrix.dat');  %各节点的连接情况
B = load('grid_center_inocean.dat');  %各节点的位置

%% 节点位置绘图
scatter(B(:,1),B(:,2));
box on;
axis equal tight;

%% 邻接矩阵的权分布
figure();
AA = A;
AA(AA == 0) = NaN;
histogram(AA);
title('degree contribution in SCS');
xlabel('weight');
ylabel('frequence');
% his_plot(1,'his_0001.nc')

%% 计算无权入度和无权出度
N = size(A);
N = N(1);
out_degree = linspace(0,0,N);
in_degree = linspace(0,0,N);

for i = 1:N
    a = find(A(i,:));
    out_degree(i) = length(a);
    b = find(A(:,i));
    in_degree(i) = length(b);
end
general_degree = out_degree + in_degree;
disp('计算完成');
A2 = A*A;


%% 计算入聚集系数和出聚集系数
% 方法一
% 1.得到每点的度
% 2.得到关联该点的所有节点的编号
% 3.在第二步得到的节点中, 那些节点之间的连接情况. 
% 使用切片取出这些节点构成的子邻接矩阵
% 使用find方法找到总的弧的数目
% 4. 实际弧数除以最大可能弧数
clustering = linspace(0,0,N); % 无向底图的聚集系数
in_clustering = linspace(0,0,N); % 入集团聚集系数
out_clustering = linspace(0,0,N); % 出集团聚集系数
for node = 1:N
    if general_degree(node) ~= 0
        in_incident = find(A(:,node));
        out_incident = find(A(node,:));
        incident = [out_incident,in_incident'];
        
        in_subA = A(in_incident',in_incident');
        out_subA = A(out_incident,out_incident);
        subA = A(incident,incident);
        
        in_arch_counts = length(find(in_subA));
        out_arch_counts = length(find(out_subA));
        arch_counts = length(find(subA));
        
        in_node_counts = in_degree(node);
        out_node_counts = out_degree(node);
        node_counts = general_degree(node); 
        
        clustering(node) = arch_counts/(node_counts*(node_counts-1));
        in_clustering(node) = in_arch_counts/(in_node_counts*(in_node_counts-1));
        out_clustering(node) = out_arch_counts/(out_node_counts*(out_node_counts-1));
    else
        clustering(node) = 0;
        in_clustering(node) = 0;
        out_clustering(node) = 0;
    end
end
length(find(clustering>1))
length(find(in_clustering>1))
length(find(out_clustering>1))


% 方法二
% A_A = A+A';
% A_A3 = A_A*A_A*A_A;
% A_A2 = A_A*A_A;
% cc = linspace(0,0,N);
% for node = 1:N
%     cc(node) = A_A3(node,node)/(A_A2(node,node)*(A_A2(node,node)-1));
% end

%% 集聚系数 绘图
figure();
d = clustering';
length(find(d))/length(d) 
d(d>1) = NaN;
scatter(B(:,1),B(:,2),69,d,'square','filled')
axis equal
box on
title('general clustering coefficient SCS');
xlabel('lon');
ylabel('lat');
colorbar

figure();
d = out_clustering';
length(find(d))/length(d) 
d(d>1) = NaN;
scatter(B(:,1),B(:,2),69,d,'square','filled')
axis equal
box on
title('out clustering coefficient SCS');
xlabel('lon');
ylabel('lat');
colorbar

figure();
d = in_clustering';
length(find(d))/length(d) 
d(d>1) = NaN;
scatter(B(:,1),B(:,2),69,d,'square','filled')
axis equal
box on
title('in clustering coefficient SCS');
xlabel('lon');
ylabel('lat');
colorbar

%% 无权度 绘图
% 入度
figure();
d = in_degree';
length(find(d))/length(d) %非零格点的占比
% d3(out_d3>100) = NaN;
scatter(B(:,1),B(:,2),69,d,'square','filled')
axis equal
box on
title('in degree SCS');
xlabel('lon');
ylabel('lat');
colorbar
% 出度
figure();
d = out_degree';
length(find(d))/length(d) %非零格点的占比
% d3(out_d3>100) = NaN;
scatter(B(:,1),B(:,2),69,d,'square','filled')
axis equal
box on
title('out degree SCS');
xlabel('lon');
ylabel('lat');
colorbar
% 总度
figure();
d = general_degree';
length(find(d))/length(d) %非零格点的占比
% d3(out_d3>100) = NaN;
scatter(B(:,1),B(:,2),69,d,'square','filled')
axis equal
box on
title('general degree SCS');
xlabel('lon');
ylabel('lat');
colorbar
figure()
% histogram(d)
% data = load('Ini2Final.dat');
% ini_grid = data(:,1);
% final_grid = data(:,2);
% histogram(ini_grid,1192)

%% 始末位置
ini = load('initial_valid.dat');                                            
final = load('final_valid.dat');

figure();
scatter(ini(1:1:end,1),ini(1:1:end,2),0.01,'black');
axis equal tight;
box on;

figure()
scatter(final(1:1:end,1),final(1:1:end,2),0.01,'black');
axis equal tight;
box on;

%%
a = load('Ini2Final.dat')
histogram(a(:,1))

%% 最可能路径绘图
pp = load('path.dat');
hold on;
plot(pp(:,1),pp(:,2));