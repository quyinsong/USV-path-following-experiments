% test 1
clc
clear all
close all

ts = 0.04;
tfinal = 720;
Ns = tfinal/ts;

for k = 1:Ns
    t = (k-1)*ts;
    % 初始化
    if t==0
        x10 = [0.1 0.1 0.1 10 0 -pi/6]'; Tc1 = [0 0]';
        x20 = [0.1 0.1 0.1 22 -1 0]'; Tc2 = [0 0]';
        x30 = [0.1 0.1 0.1 29 -1.5 pi]'; Tc3 = [0 0]';
        Tmax = 113.6775; Tmin = -55;
        Vc = [0 0]';
        d11 = 151.57;
        hT1 = zeros(2,1); hT2 = zeros(2,1); hT3 = zeros(2,1);
    end
    % 外界干扰
    tau_w = [5*sin(0.08*t)*cos(0.15*t)+5 3*cos(0.08*t)*sin(0.15*t)+3 ...
            5*cos(0.08*t)*sin(0.15*t)+5]';
    % 无人艇模型
    [x1, xdot1, T1, fu1, fr1] = ASV1( x10, Tc1, tau_w, Vc, ts);
    [x2, xdot2, T2, fu2, fr2] = ASV2( x20, Tc2, tau_w, Vc, ts);
    [x3, xdot3, T3, fu3, fr3] = ASV3( x30, Tc3, tau_w, Vc, ts);
    
    % 传感器获得值，及相应计算量
    beta1 = atan2(x1(2),x1(1));  % 漂角
    beta2 = atan2(x2(2),x2(1));  % 漂角
    beta3 = atan2(x3(2),x3(1));  % 漂角
    psi1 = x1(6); psi2 = x2(6); psi3 = x3(6);% 艏向角
    u1 = x1(1); v1 = x1(2); r1 = x1(3); % 无人艇速度和角速度状态
    u2 = x2(1); v2 = x2(2); r2 = x2(3); % 无人艇速度和角速度状态
    u3 = x3(1); v3 = x3(2); r3 = x3(3); % 无人艇速度和角速度状态
    U1 = sqrt(u1^2+v1^2);      % 无人艇总航向速度
    U2 = sqrt(u2^2+v2^2);      % 无人艇总航向速度
    U3 = sqrt(u3^2+v3^2);      % 无人艇总航向速度
    
    % communication networks
    if t==0
        A = [0 0 0;
             1 0 0;
             0 1 0;];
        B = diag([1 0 0]);
        L = [0 0 0;
            -1 1 0;
             0 -1 1;];
        Dw = [0 0 0;
              0 0 0;
              0 0 0;];
        Dw0 = [0;0;0];
        w = [0.1 0.8 0.5]';
        us0 = 0.3;
        w0 = 1.5;
    end
    
    % 混合路径
    xd0 = 10;
    yd0 = w0;
    xd0_dw0 = 0;
    yd0_dw0 = 1;
    if w0 >= 80
         R = 10;
         xd0 = R*cos(w0-80);
         yd0 = R*sin(w0-80)+80;
         xd0_dw0 = -R*sin(w0-80);
         yd0_dw0 = R*cos(w0-80);
     end
     if w0 >= 80+pi
         xd0 = -10;
         yd0 = -(w0-80-pi)+80;
         xd0_dw0 = 0;
         yd0_dw0 = -1;
     end
     if w0 >= 160+pi
         xd0 = R*cos(w0-160-pi+pi);
         yd0 = R*sin(w0-160-pi+pi);
         xd0_dw0 = -R*sin(w0-160-pi+pi);
         yd0_dw0 = R*cos(w0-160-pi+pi);
     end
     if w0 >= 160+2*pi
         xd0 = 10;
         yd0 = w0-160-2*pi;
         xd0_dw0 = 0;
         yd0_dw0 = 1;
     end
     
     % 圆形路径
%      R = 10;
%      xd0 = R*cos(w0); 
%      yd0 = R*sin(w0); 
     
     % 正弦路径
     
    vs = us0/sqrt(xd0_dw0^2+yd0_dw0^2);  
    w0_dot = vs;
    w0 = ts*w0_dot+w0;
    
    ew = (L+B)*w-sum(A.*Dw,2)-B*(Dw0+w0.*ones(3,1));
    % LOS制导
    eta1 = [x1(4),x1(5),x1(6)]'; eta2 = [x2(4),x2(5),x2(6)]'; eta3 = [x3(4),x3(5),x3(6)]';
    nu1 = [x1(1),x1(2),x1(3)]'; nu2 = [x2(1),x2(2),x2(3)]';  nu3 = [x3(1),x3(2),x3(3)]'; 
    [xd1, yd1, alpha_psi1, psip1, alpha_r1, alpha_u1, xe1, ye1, w(1), a1, b1, psie1] = LOSguidance1( 4, 5, eta1, nu1, vs, w(1), ts, t, ew(1) );
    [xd2, yd2, alpha_psi2, psip2, alpha_r2, alpha_u2, xe2, ye2, w(2), a2, b2, psie2] = LOSguidance2( 4, 5, eta2, nu2, vs, w(2), ts, t, ew(2) );
    [xd3, yd3, alpha_psi3, psip3, alpha_r3, alpha_u3, xe3, ye3, w(3), a3, b3, psie3] = LOSguidance3( 4, 5, eta3, nu3, vs, w(3), ts, t, ew(3) );
    
    
   %利用双曲正弦函数近似饱和非线性
   for i=1:2
       if(Tc1(i) > 0)
           hT1(i) = Tmax*tanh(Tc1(i)/Tmax); 
       end
      if(Tc1(i)<=0)
           hT1(i) = Tmin*tanh(Tc1(i)/Tmin); 
      end
      if(Tc2(i) > 0)
           hT2(i) = Tmax*tanh(Tc2(i)/Tmax); 
      end
      if(Tc2(i)<=0)
           hT2(i) = Tmin*tanh(Tc2(i)/Tmin); 
      end
      if(Tc3(i) > 0)
           hT3(i) = Tmax*tanh(Tc3(i)/Tmax); 
      end
      if(Tc3(i)<=0)
           hT3(i) = Tmin*tanh(Tc3(i)/Tmin); 
      end  
   end
   
   % NN controller
   wd1 = [alpha_u1,alpha_r1]'; wd2 = [alpha_u2,alpha_r2]'; wd3 = [alpha_u3,alpha_r3]';
   w1 = [u1,r1]'; w2 = [u2,r2]'; w3 = [u3,r3]';
   [Fhat1, lambda1, Tc1, ew1] = NNcontroller1( wd1,w1,hT1,Tc1,v1,ts );
   [Fhat2, lambda2, Tc2, ew2] = NNcontroller2( wd2,w2,hT2,Tc2,v2,ts );
   [Fhat3, lambda3, Tc3, ew3] = NNcontroller3( wd3,w3,hT3,Tc3,v3,ts );
   
   % 存储时间序列
   tout(k,1) = t;
   Tlim(k,:) = [Tmax,Tmin];
   xout1(k,:) = [x1',Tc1',T1',u1];
   xout2(k,:) = [x2',Tc2',T2',u2];
   xout3(k,:) = [x3',Tc3',T3',u3];
   Eout(k,:) = [ew',xe1,xe2,xe3,ye1,ye2,ye3];
   eout(k,:) = [psie1, psie2, psie3, ew1', ew2', ew3'];
   Pdout(k,:) = [xd1,yd1,xd2,yd2,xd3,yd3,xd0,yd0]';
   Fout(k,:) = [fu1,fu2,fu3,fr1,fr2,fr3,Fhat1(1),Fhat2(1),Fhat3(1),Fhat1(2),Fhat2(2),Fhat3(2)]';
   psipout(k,:) = [psip1,psip2,psip3];
   wout(k,:) = [w0,w'];
   about(k,:) = [a1,b1,a2,b2,a3,b3];
   
   hT(k,:) = [hT1',hT2',hT3']';
end
t = tout(:,1);

u1 = xout1(:,1); v1 = xout1(:,2); r1 = xout1(:,3); x1 = xout1(:,4);
y1 = xout1(:,5); psi1 = xout1(:,6); Tc1r = xout1(:,7); Tc1l = xout1(:,8);
T1r = xout1(:,9); T1l = xout1(:,10); u1 = xout1(:,11);

u2 = xout2(:,1); v2 = xout2(:,2); r2 = xout2(:,3); x2 = xout2(:,4);
y2 = xout2(:,5); psi2 = xout2(:,6); Tc2r = xout2(:,7); Tc2l = xout2(:,8);
T2r = xout2(:,9); T2l = xout2(:,10); u2 = xout2(:,11);

u3 = xout3(:,1); v3 = xout3(:,2); r3 = xout3(:,3); x3 = xout3(:,4);
y3 = xout3(:,5); psi3 = xout3(:,6); Tc3r = xout3(:,7); Tc3l = xout3(:,8);
T3r = xout3(:,9); T3l = xout3(:,10); u3 = xout3(:,11);

% plots
linewid = 1;
fontsize = 12;
fontname = 'Time news roman';
fontweight = 'bold';

%绘制协同路径跟踪效果
figure(1);hold on
box on

xrange=[-50 20 120]; yrange = [-50 10 50];
for k=1:1000:Ns
    pos1 = [x1(k),y1(k)]'; pos2 = [x2(k),y2(k)]'; pos3 = [x3(k),y3(k)]';
    plotCircle( 2.5, [Pdout(k,8),Pdout(k,7)], 'k-', linewid );
    modelplot(pos1,psi1(k),xrange,yrange,'r-',linewid);
    modelplot(pos2,psi2(k),xrange,yrange,'b-',linewid);
    modelplot(pos3,psi3(k),xrange,yrange,'g-',linewid);
end
Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);
plt1 = plot(y1,x1,'r-','linewidth',linewid);
plt2 = plot(y2,x2,'b-','linewidth',linewid);
plt3 = plot(y3,x3,'g-','linewidth',linewid);

plt4 = plot(Pdout(:,2),Pdout(:,1),'r--','linewidth',2);
plt5 = plot(Pdout(:,4),Pdout(:,3),'b--','linewidth',2);
plt6 = plot(Pdout(:,6),Pdout(:,5),'g--','linewidth',2);

plt7 = plot(Pdout(:,8),Pdout(:,7),'k--','linewidth',2);

h1 = xlabel('$$y\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('$$x\ $$ (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

% h1=legend([plt1(1),plt2(1),plt3(1)],{'trajectory1','trajectory2','trajectory3'});
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% set(h1,'FontWeight',fontweight); set(h1,'Box','off'); set(h1,'orientation','horizontal');
% 
% ah=axes('position',get(gca,'position'),'visible','off');
% h2 = legend(ah,[plt4(1),plt5(1),plt6(1)],'desired path1','desired path2','desired path3');
% set(h2,'Interpreter','latex'); set(h2,'FontSize',fontsize); set(h2,'FontName',fontname);
% set(h2,'FontWeight',fontweight); set(h2,'Box','off'); set(h2,'orientation','horizontal');
% 
% ah=axes('position',get(gca,'position'),'visible','off');
% h3 = legend(ah,plt7(1),'virtual path');
% set(h3,'Interpreter','latex'); set(h3,'FontSize',fontsize); set(h3,'FontName',fontname);
% set(h3,'FontWeight',fontweight); set(h3,'Box','off'); set(h3,'orientation','horizontal');
% box off;

h1=legend([plt1(1),plt4(1),plt7(1)],{'trajectory','desired path','virtual path'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal');

set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);
hold off;


tmin = 0;
tmax = tfinal;
tinterval = tfinal/6;

%绘制误差约束效果
figure(2)
subplot(211)
plot(t,Eout(:,4),'r-',t,Eout(:,5),'b-',t,Eout(:,6),'g-','linewidth',linewid); hold on
plot(t,about(:,1),'c--',t,about(:,3),'m--',t,about(:,5),'k--','linewidth',linewid);
plot(t,-about(:,1),'c--',t,-about(:,3),'m--',t,-about(:,5),'k--','linewidth',linewid); hold off
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
h1=legend('USV1','USV2','USV3');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);  set(h1,'orientation','horizontal');
xlabel('time (s)');
h1 = ylabel('$$x_{e}\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
axes('Position',[0.2,0.85,0.2,0.1]); % 生成子图   
plot(t,Eout(:,4),'r-',t,Eout(:,5),'b-',t,Eout(:,6),'g-','linewidth',linewid); hold on
plot(t,about(:,1),'c--',t,about(:,3),'m--',t,about(:,5),'k--','linewidth',linewid);
plot(t,-about(:,1),'c--',t,-about(:,3),'m--',t,-about(:,5),'k--','linewidth',linewid); hold off; % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([-5,5]);
set(gca,'xTick',0:5:20);
set(gca,'yTick',-5:5:5);


subplot(212)
plot(t,Eout(:,7),'r-',t,Eout(:,8),'b-',t,Eout(:,9),'g-','linewidth',linewid); hold on
plot(t,about(:,2),'c--',t,about(:,4),'m--',t,about(:,6),'k--','linewidth',linewid);
plot(t,-about(:,2),'c--',t,-about(:,4),'m--',t,-about(:,6),'k--','linewidth',linewid); hold off
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
xlabel('time (s)');
h1 = ylabel('$$y_{e}\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
axes('Position',[0.2,0.38,0.2,0.1]); % 生成子图   
plot(t,Eout(:,7),'r-',t,Eout(:,8),'b-',t,Eout(:,9),'g-','linewidth',linewid); hold on
plot(t,about(:,2),'c--',t,about(:,4),'m--',t,about(:,6),'k--','linewidth',linewid);
plot(t,-about(:,2),'c--',t,-about(:,4),'m--',t,-about(:,6),'k--','linewidth',linewid); hold off; % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([-5,5]);
set(gca,'xTick',0:5:20);
set(gca,'yTick',-5:5:5);

%绘制干扰估计效果
figure(3)
subplot(211)
plot(t,Fout(:,1),'r-',t,Fout(:,7),'b--','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
h1 = legend('$$f_{u1}$$','$$\hat{f}_{u1}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');
h1 = ylabel('$$f_{u1}\ $$ (N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
axes('Position',[0.2,0.68,0.2,0.15]); % 生成子图   
plot(t,Fout(:,1),'r-',t,Fout(:,7),'b--','linewidth',linewid); % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([-5,0]);
set(gca,'xTick',0:5:20);
set(gca,'yTick',-5:2.5:0);

subplot(212)
plot(t,Fout(:,4),'r-',t,Fout(:,10),'b--','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
h1 = legend('$$f_{r1}$$','$$\hat{f}_{r1}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');
h1 = ylabel('$$f_{r1}\ $$ (Nm)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
axes('Position',[0.2,0.18,0.2,0.15]); % 生成子图   
plot(t,Fout(:,4),'r-',t,Fout(:,10),'b--','linewidth',linewid); % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([-8,2]);
set(gca,'xTick',0:5:20);
set(gca,'yTick',-8:5:2);

%绘制三个速度变化
figure(4)
subplot(311)
plot(t,u1,'r-',t,u2,'b-',t,u3,'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([0 1.5]);
h1 = legend('USV1','USV2','USV3');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');
h1 = ylabel('$$u\ $$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(312)
plot(t,v1,'r-',t,v2,'b-',t,v3,'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([-0.2 0.2]);
xlabel('time (s)');
h1 = ylabel('$$v\ $$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(313)
plot(t,r1,'r-',t,r2,'b-',t,r3,'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([-0.8 0.8]);
xlabel('time (s)');
h1 = ylabel('$$r\ $$ (rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

%绘制航向角变化
figure(5)
plot(t,psipout(:,1),'r-',t,psipout(:,2),'b-',t,psipout(:,3),'g--','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);

%绘制推力饱和处理效果
figure(6)
plot(t,Tc1r,'r-',t,Tc1l,'b-','linewidth',linewid); hold on
plot(t,T1r,'g--',t,T1l,'m--','linewidth',linewid);
plot(t,Tlim(:,1),'c-',t,Tlim(:,2),'k--','linewidth',linewid); hold off
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([-350 850]);
xlabel('time (s)');
ylabel('Thrust (N)');
h1 = legend('$$T_{c1}$$','$$T_{c2}$$','$$T_{1}$$','$$T_{2}$$','$$T_{max}$$','$$T_{min}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal');
axes('Position',[0.2,0.6,0.3,0.3]); % 生成子图   
plot(t,Tc1r,'r-',t,Tc1l,'b-','linewidth',linewid); hold on
plot(t,T1r,'g--',t,T1l,'m--','linewidth',linewid);
plot(t,Tlim(:,1),'c-',t,Tlim(:,2),'k--','linewidth',linewid); hold off % 绘制局部曲线图                        
xlim([0,50]); % 设置坐标轴范围  
ylim([-120,180]);
set(gca,'xTick',0:10:50);
set(gca,'yTick',-120:60:180);

axes('Position',[0.6,0.6,0.3,0.3]); % 生成子图   
plot(t,Tc1r,'r-',t,Tc1l,'b-','linewidth',linewid); hold on
plot(t,T1r,'g--',t,T1l,'m--','linewidth',linewid);
plot(t,Tlim(:,1),'c-',t,Tlim(:,2),'k--','linewidth',linewid); hold off % 绘制局部曲线图                        
xlim([260,280]); % 设置坐标轴范围  
ylim([-120,180]);
set(gca,'xTick',260:10:280);
set(gca,'yTick',-120:60:180);

%绘制推力变化曲线
figure(7)
subplot(311)
plot(t,T1r,'r-',t,T1l,'b-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
subplot(312)
plot(t,T2r,'r-',t,T2l,'b-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
subplot(313)
plot(t,T3r,'r-',t,T3l,'b-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);

%绘制随时间变化的路径参数
figure(8)
plot(t,wout(:,1),'r-',t,wout(:,2),'g-',t,wout(:,3),'b-',t,wout(:,4),'c-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
xlabel('time (s)');
h1 = ylabel('$$\theta\ $$ (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
axes('Position',[0.2,0.6,0.25,0.25]); % 生成子图   
plot(t,wout(:,1),'r-',t,wout(:,2),'g-',t,wout(:,3),'b-',t,wout(:,4),'c-','linewidth',linewid) % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([0,10]);
set(gca,'xTick',0:5:20);
set(gca,'yTick',0:5:10);

%绘制虚拟信号跟踪误差
figure(9)
subplot(311)
plot(t,eout(:,1),'r-',t,eout(:,2),'b-',t,eout(:,3),'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([-1.5 1.5]);
h1 = legend('ASV1','ASV2','ASV3');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);  set(h1,'orientation','horizontal');
h1 = ylabel('$$\psi_e\ $$ (rad)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(312)
plot(t,eout(:,5),'r-',t,eout(:,7),'b-',t,eout(:,9),'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([-5 5]);
h1 = ylabel('$$r_e\ $$ (rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

subplot(313)
plot(t,eout(:,4),'r-',t,eout(:,6),'b-',t,eout(:,8),'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
ylim([-4 2]);
h1 = ylabel('$$u_e\ $$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');

%绘制编队误差
figure(10)
plot(t,Eout(:,1),'r-',t,Eout(:,2),'b-',t,Eout(:,3),'g-','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);
h1 = legend('USV1','USV2','USV3');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); set(h1,'orientation','horizontal');
xlabel('time (s)');
h1 = ylabel('$$e_{\theta}\ $$ (rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
axes('Position',[0.3,0.3,0.4,0.2]); % 生成子图   
plot(t,Eout(:,1),'r-',t,Eout(:,2),'b-',t,Eout(:,3),'g-','linewidth',linewid) % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([-2,1]);
set(gca,'xTick',0:5:20);
set(gca,'yTick',-2:1:1);

%绘制tanh饱和近似效果
figure(11)
plot(t,T1r,'r-',t,hT(:,1),'b-',t,T1l,'g-',t,hT(:,2),'m','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);

h1 = legend('T1r','hT1r','T1l','hT1l');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1 = ylabel('$$u_e\ $$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');

figure(12)
plot(t,T2r,'r-',t,hT(:,3),'b-',t,T2l,'g-',t,hT(:,4),'m','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);

h1 = legend('T2r','hT2r','T2l','hT2l');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1 = ylabel('$$u_e\ $$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');

figure(13)
plot(t,T3r,'r-',t,hT(:,5),'b-',t,T3l,'g-',t,hT(:,6),'m','linewidth',linewid)
set(gca,'xTick',tmin:tinterval:tmax);
xlim([tmin tmax]);

h1 = legend('T3r','hT3r','T3l','hT3l');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

h1 = ylabel('$$u_e\ $$ (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlabel('time (s)');














