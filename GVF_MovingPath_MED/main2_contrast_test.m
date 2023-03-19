clc
clear all
close all

load data4
xout1 = xout;
load data5
xout2 = xout;
load data6
xout3 = xout;
clear xout

y1 = xout1(:,4); x1 = xout1(:,5); psi1 = xout1(:,6);
t1 = xout1(:,7);
yT1 = xout1(:,8); xT1 = xout1(:,9);
phi1 = xout1(:,10); wd1 = xout1(:,11); 
w1 = xout1(:,16); what1 = xout1(:,17); 
alphaT1 = xout1(:,27);

y2 = xout2(:,4); x2 = xout2(:,5); psi2 = xout2(:,6);
yT2 = xout2(:,8); xT2 = xout2(:,9);
phi2 = xout2(:,10); wd2 = xout2(:,11); 
w2 = xout2(:,16); what2 = xout2(:,17); 
alphaT2 = xout2(:,27);

y3 = xout3(:,4); x3 = xout3(:,5); psi3 = xout3(:,6);
t3 = xout3(:,7);
yT3 = xout3(:,8); xT3 = xout3(:,9);
phi3 = xout3(:,10); wd3 = xout3(:,11); 
w3 = xout3(:,16); what3 = xout3(:,17); 
alphaT3 = xout3(:,27);
%% PLOTS
ts = 0.02;
tfinal = 120;
Ns = tfinal/ts;

fontsize = 14; fontname = 'TimesNewRoman'; linewid = 2; fontweight = 'bold';

figure(1); hold on
xrange=[-25 10 25]; yrange = [-5 10 55];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

xd = 0:0.2:50;
yd = 20*sin(2*pi*xd/50);
plot(yd,xd,'k-','linewid',2);
plot(x1,y1,'r--','linewid',2);
plot(x2,y2,'g--','linewid',2);
plot(x3,y3,'b--','linewid',2);

% plot(xT,yT,'b-','linewid',2);
h1=legend('desired path','$$trajectory\ k_{mp}=0.05$$',...
    '$$trajectory\ k_{mp}=0.1$$','$$trajectory\ k_{mp}=0.4$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% y_tick = 13.5784; x_tick = 13.5521;
% text(x_tick+1.5,y_tick,'\leftarrow Faults occur','color','r','FontSize',12,'linewid',2.5);
% plot(x_tick,y_tick,'*','MarkerSize',15,'color','r')

for k=1:800:Ns
    pos =[y1(k) x1(k)]';
    modelplot(pos,psi1(k),xrange,yrange);
end

for k=1:800:Ns
    pos =[y2(k) x2(k)]';
    modelplot1(pos,psi2(k),xrange,yrange);
end

for k=1:800:Ns
    pos =[y3(k) x3(k)]';
    modelplot2(pos,psi3(k),xrange,yrange);
end

% set(gca,'xTick',Xmin:Xinterval:Xmax);
% set(gca,'yTick',Ymin:Yinterval:Ymax);
% axis([Xmin Xmax,Ymin Ymax]);

h1=xlabel('$$y (m)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$x (m)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

% for i=1:NO
% %     xoi = O(i,1); yoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
%     yoi = O(i,1); xoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
%     h=rectangle('Position',[xoi-Roi,yoi-Roi,2*Roi,2*Roi],'Curvature',[1,1],'EdgeColor','k');
%     set(h,'LineStyle','--','linewid',2);
%     h=rectangle('Position',[xoi-Ro1i,yoi-Ro1i,2*Ro1i,2*Ro1i],'Curvature',[1,1],'EdgeColor','k');
%     set(h,'LineStyle','--','linewid',2);
% end

hold off

figure(2); hold on
plot(t1,phi1,'r-',t1,phi2,'g-',t3,phi3,'b-','linewid',2);
h1=ylabel('$$\phi(\xi) (m^{2})$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = legend('$$k_{mp}=0.05$$','$$k_{mp}=0.1$$','$$k_{mp}=0.4$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
hold off

figure(3); 
subplot(3,1,1);plot(t1,wd1,'r-',t1,w1,'b-',t1,what1,'c--','linewid',2);
h1 = legend('$$w_{d}$$','$$w$$','$$\hat{w}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% h1 = xlabel('t (s)');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% h1 = ylabel('(rad/s)');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
h1 = title('$$k_{mp}=0.05$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
subplot(3,1,2);plot(t1,wd2,'r-',t1,w2,'b-',t1,what2,'c--','linewid',2);
% h1 = legend('$$w_{d}$$','$$w$$','$$\hat{w}$$');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% h1 = xlabel('t (s)');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
h1 = ylabel('(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
h1 = title('$$k_{mp}=-0.1$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
subplot(3,1,3);plot(t1,wd3,'r-',t1,w3,'b-',t1,what3,'c--','linewid',2);
% h1 = legend('$$w_{d}$$','$$w$$','$$\hat{w}$$');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
h1 = xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
h1 = title('$$k_{mp}=0.4$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% h1 = ylabel('(rad/s)');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);




