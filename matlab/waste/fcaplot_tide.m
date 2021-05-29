%% 入度
var = s1;
scatter(var(:,1),var(:,2))
xticks(linspace(-0.005,0.03,8))
yticks(linspace(0,30,7))
xlabel('backawrd FTLE');
ylabel('in degree');



axis square
    axis tight
hold on
plot([-0.005,0.03],[0,30],'r')

box on

figname = 'sheet2.png'
fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)

%% m2a
ddda = m2a
xcontent = '实测值/cm'
ycontent = '计算值/cm'
figname = 'M2振幅.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)
xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 140

xticks(linspace(0,140,8))
xticklabels(linspace(0,140,8))
yticks(linspace(0,140,8))
yticklabels(linspace(0,140,8))
axis equal

xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 4.72cm','FontSize',14)
text(bigger*0.1,bigger*0.92,'M_2振幅','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)

%% s2a
ddda = s2a
xcontent = '实测值/cm'
ycontent = '计算值/cm'
figname = 'S2振幅.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)
xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 40

xticks(linspace(0,40,5))
xticklabels(linspace(0,40,5))
yticks(linspace(0,40,5))
yticklabels(linspace(0,40,5))
axis equal

xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 1.42cm','FontSize',14)
text(bigger*0.1,bigger*0.92,'S_2振幅','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)

%% k1a
ddda = k1a
xcontent = '实测值/cm'
ycontent = '计算值/cm'
figname = 'K1振幅.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)
xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 50

xticks(linspace(0,50,6))
xticklabels(linspace(0,50,6))
yticks(linspace(0,50,6))
yticklabels(linspace(0,50,6))
axis equal

xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 3.75cm','FontSize',14)
text(bigger*0.1,bigger*0.92,'K_1振幅','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)

%% o1a
ddda = o1a
xcontent = '实测值/cm'
ycontent = '计算值/cm'
figname = 'O1振幅.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)
xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 30

xticks(linspace(0,30,4))
xticklabels(linspace(0,30,4))
yticks(linspace(0,30,4))
yticklabels(linspace(0,30,4))
axis equal

xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 2.95cm','FontSize',14)
text(bigger*0.1,bigger*0.92,'0_1振幅','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)







%% m2p
ddda = m2p
xcontent = '实测值/°'
ycontent = '计算值/°'
figname = 'M2迟角.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)

xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 360

xticks(linspace(0,360,7))
xticklabels(linspace(0,360,7))
yticks(linspace(0,360,7))
yticklabels(linspace(0,360,7))

axis equal
xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 10.06°','FontSize',14)
text(bigger*0.1,bigger*0.92,'M_2迟角','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)




%% s2p
ddda = s2p
xcontent = '实测值/°'
ycontent = '计算值/°'
figname = 'S2迟角.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)

xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 360

xticks(linspace(0,360,7))
xticklabels(linspace(0,360,7))
yticks(linspace(0,360,7))
yticklabels(linspace(0,360,7))

axis equal
xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 11.46°','FontSize',14)
text(bigger*0.1,bigger*0.92,'S_2迟角','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)




%% k1p
ddda = k1p
xcontent = '实测值/°'
ycontent = '计算值/°'
figname = 'K1迟角.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)

xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 360

xticks(linspace(0,360,7))
xticklabels(linspace(0,360,7))
yticks(linspace(0,360,7))
yticklabels(linspace(0,360,7))

axis equal
xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 9.42°','FontSize',14)
text(bigger*0.1,bigger*0.92,'K_1迟角','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)



%% o1p
ddda = o1p
xcontent = '实测值/°'
ycontent = '计算值/°'
figname = 'O1迟角.png'

scatter(ddda(:,1),ddda(:,2),'filled','black')
set(gca,'fontsize',14)

xmax = max(ddda(:,1))
ymax = max(ddda(:,2))
bigger = max(xmax,ymax)
bigger = 360

xticks(linspace(0,360,7))
xticklabels(linspace(0,360,7))
yticks(linspace(0,360,7))
yticklabels(linspace(0,360,7))

axis equal
xlabel(xcontent)
ylabel(ycontent)
xlim([0,bigger])
ylim([0,bigger])
box on
text(bigger*0.1,bigger*0.85,'MAE = 5.01°','FontSize',14)
text(bigger*0.1,bigger*0.92,'0_1迟角','FontSize',14)
%title(titlecontent,'FontSize',14)
hold on 
plot(0:bigger,0:bigger,'--')

fig = gcf;     
print(fig, '-dpng', '-r600',figname)   
hold off
clf    
close(fig)

