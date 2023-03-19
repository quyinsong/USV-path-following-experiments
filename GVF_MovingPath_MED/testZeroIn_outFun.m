clc
clear all
close all

c0 = -(1-20)^2;% c = -361
x0=(c0+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
d1 = -(c0-x0); % 距离障碍物中心距离，越小距离障碍物越近
d2 = -x0; % 距离反应域边界距离，越大距离障碍物越远
l1 = 10;
l2 = 0.1;
fx0 = exp(l1./(c0-x0)); % 路径跟踪权值
fx02 = exp(l2./x0); % 避障权值

c1 = -(5-20)^2;% c = -255
x1=(c1+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
fx1 = exp(l1./(c1-x1)); % 路径跟踪权值
fx12 = exp(l2./x1); % 避障权值

c2 = -(10-20)^2;% c = -100
x2=(c2+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
fx2 = exp(l1./(c2-x2)); % 路径跟踪权值
fx22 = exp(l2./x2); % 避障权值

c3 = -(15-20)^2;% c = -25
x3=(c3+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
fx3 = exp(l1./(c3-x3)); % 路径跟踪权值
fx32 = exp(l2./x3); % 避障权值

c4 = -(19-20)^2;% c = -1
x4=(c4+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
fx4 = exp(l1./(c4-x4)); % 路径跟踪权值
fx42 = exp(l2./x4); % 避障权值

%% PLOTS
fontsize = 14; fontname = 'TimesNewRoman'; linewid = 2; fontweight = 'bold';

figure(1);hold on
% plot(x,fx,'r--',x,fx1,'b--','linewid',2);
plot(x0,fx0./(fx0+fx02),'r-',x0,fx02./(fx0+fx02),'r--','linewid',2);
plot(x1,fx1./(fx1+fx12),'g-',x1,fx12./(fx1+fx12),'g--','linewid',2);
plot(x2,fx2./(fx2+fx22),'b-',x2,fx22./(fx2+fx22),'b--','linewid',2);
plot(x3,fx3./(fx3+fx32),'k-',x3,fx32./(fx3+fx32),'k--','linewid',2);
plot(x4,fx4./(fx4+fx42),'m-',x4,fx42./(fx4+fx42),'m--','linewid',2);
legend('路径跟踪权值 Q=1','避障权值 Q=1','路径跟踪权值 Q=5','避障权值 Q=5','路径跟踪权值 Q=10','避障权值 Q=10',...
       '路径跟踪权值 Q=15','避障权值 Q=15','路径跟踪权值 Q=19','避障权值 Q=19');
line([-399,-399],[0,1],'color',[1 1 0]); 
line([-336,-336],[0,1],'color',[1 0 0]);
line([-300,-300],[0,1],'color',[0 1 0]);
line([-175,-175],[0,1],'color',[0 0 1]);
line([-39,-39],[0,1],'color',[0 1 1]);
hold off
   
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------

c0 = -(10-20)^2;% c = -100
x0=(c0+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
d1 = -(c0-x0); % 距离障碍物中心距离，越小距离障碍物越近
d2 = -x0; % 距离反应域边界距离，越大距离障碍物越远
l0 = 0.1;
l02 = 1;
fx0 = exp(l0./(c0-x0)); % 路径跟踪权值
fx02 = exp(l02./x0); % 避障权值

c1 = -(10-20)^2;% c = -100
x1=(c1+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
l1 = 1;
l12 = 1;
fx1 = exp(l1./(c1-x1)); % 路径跟踪权值
fx12 = exp(l12./x1); % 避障权值

c2 = -(10-20)^2;% % c = -100
x2=(c2+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
l2 = 10;
l22 = 1;
fx2 = exp(l2./(c2-x2)); % 路径跟踪权值
fx22 = exp(l22./x2); % 避障权值

c3 = -(10-20)^2;% c = -100
x3=(c3+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
l3 = 100;
l32 = 1;
fx3 = exp(l3./(c3-x3)); % 路径跟踪权值
fx32 = exp(l32./x3); % 避障权值

c4 = -(10-20)^2; % c = -100;
x4=(c4+0.00001):0.01:-0.00001; % when x=-336, x locate in re of circle
l4 = 1000;
l42 = 1;
fx4 = exp(l4./(c4-x4)); % 路径跟踪权值
fx42 = exp(l42./x4); % 避障权值

figure(2);hold on
% plot(x,fx,'r--',x,fx1,'b--','linewid',2);
plot(x1,fx1./(fx1+fx12),'r-',x1,fx12./(fx1+fx12),'r--','linewid',2);
plot(x2,fx2./(fx2+fx22),'b-',x2,fx22./(fx2+fx22),'b--','linewid',2);
plot(x3,fx3./(fx3+fx32),'k-',x3,fx32./(fx3+fx32),'k--','linewid',2);
plot(x4,fx4./(fx4+fx42),'m-',x4,fx42./(fx4+fx42),'m--','linewid',2);
h1=legend('$$S_{p}: l_{p} = 1$$','$$S_{o}: l_{p} = 1$$',...
     '$$S_{p}: l_{p} = 10$$','$$S_{o}: l_{p} = 10$$',...
     '$$S_{p}: l_{p} = 100$$','$$S_{o}: l_{p} = 100$$',...
     '$$S_{p}: l_{p} = 1000$$','$$S_{o}: l_{p} = 1000$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
set(h1,'Location','EastOutside');
h1 = xlabel('$$\varphi(\xi)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$S_{p}\ and\ S_{o}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([-100 0]);
ylim([0 1]);
% line([-100,-100],[0,1],'color',[0 0.5 0.5],'linestyle','-','linewid',2);
box on
hold off
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
