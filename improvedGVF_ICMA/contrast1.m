% date: 9 29 2022
% Author: Quyinsong

% 仿真实验1 基于观测器的向量场制导实验
clc
clear all
close all

ts = 0.02;
Tfinal = 250;
Ns = Tfinal/ts;
% 两次仿真实验相同参数
pc = 2;  % 路径相同
ud = 0.6; % 纵荡速度相同
V = 0.4;  % 流速相同
w = 0.12;  % 流的方向角速度（有旋流）
for k=1:Ns
   t = (k-1)*ts;
   % 初始化
   if t==0
       psi = 90*pi/180; x = [0 0 0 10 10 psi]'; tau_c = [0 0]';
       dp = 0.26; B = [1 1;dp -dp];
       Vc = [0 0]'; 
       E = [0 -1;1 0]; I2 = diag([1 1]);
       Tmax = 113.6775; Tmin = -55;
   end
   % 海流变化
%    if t>=50
%        Vc = [0.1 0.6]';
%    end
   if t>=50
       Vc = [V*sin(w*t) V*cos(w*t)]';
   end
   % 外界干扰
   tau_w = [10*sin(0.03*t)*cos(0.2*t)+10 15*cos(0.03*t)*sin(0.2*t)...
            5*cos(0.03*t)*sin(0.2*t)+3]';
   % 无人艇模型
   Tc = B\tau_c;
   [x, xdot, T, fu, fr] = USV2( x, Tc, tau_w, Vc, ts);
   % 计算饱和差值
   tau = B*T;  % 真实值
   Delta_tau = tau_c-tau;  % 饱和差值
   
   % 传感器获得值，及相应计算量
   beta = atan2(x(2),x(1));  % 漂角
   X = x(6)+beta;            % 航向角
   psi = x(6);               % 艏向角
   p = [x(4),x(5)]';             % 无人艇位置向量
   u = x(1); v = x(2); r = x(3); % 无人艇速度和角速度状态
   U = sqrt(x(1)^2+x(2)^2);      % 无人艇总航向速度
   m1 = [cos(psi) sin(psi)]';    % 无人艇姿态向量1
   p_dot = [xdot(4), xdot(5)]';  % 无人艇位置导数向量
   m2 = [-sin(psi) cos(psi)]';   % 无人艇姿态向量2
   vc = v*m2+Vc;                 % 待估计的未知量（即侧滑速度和海流速度）
   % 海流观测器
   [ vchat,vchat_dot,phat] = observer2( p,u,psi,ts );
   
   % 向量场制导
   [rd,e,xd,yd,theta,e_alphap,Xd_dot]  = GVFguidance2( pc,ud,vchat,vchat_dot,p,p_dot,psi,ts );

   % surge controller
   [uhat, fuhat ,tuc, lambdau] = ESOcontroller1( ud,u,tau(1),Delta_tau(1), ts );
   
   % heading controller
   [rhat, frhat, trc, lambdar] = ESOcontroller2( rd,r,tau(2),Delta_tau(2),ts );
   
   % total control signal
   tau_c = [tuc,trc]';

   % 存储时间序列
   xout(k,:) = [x' t e' Vc' vchat' fu fr fuhat frhat tau' tau_c'...
                xd yd ud  rd theta Tc' T' e_alphap Xd_dot vchat_dot' phat' vc'...
                lambdau lambdar];
end

test1.u = xout(:,1); test1.v = xout(:,2); test1.r = xout(:,3);
test1.x = xout(:,4); test1.y = xout(:,5); test1.psi = xout(:,6);
test1.t = xout(:,7); test1.phi1 = xout(:,8); test1.phi2 = xout(:,9);
test1.Vc = xout(:,10:11); test1.vchat = xout(:,12:13);
test1.fu = xout(:,14); test1.fr = xout(:,15);
test1.fuhat = xout(:,16); test1.frhat = xout(:,17);
test1.tau = xout(:,18:19); test1.tau_c = xout(:,20:21);
test1.xd = xout(:,22); test1.yd = xout(:,23);
test1.ud = xout(:,24);
test1.rd = xout(:,25);
test1.theta = xout(:,26);
test1.Tc = xout(:,27:28); test1.T = xout(:,29:30);
test1.e_alphap = xout(:,31); test1.Xd_dot = xout(:,32);
test1.vchat_dot= xout(:,33:34);
test1.xhat = xout(:,35); test1.yhat = xout(:,36);
test1.vc = xout(:,37:38); 
test1.lambdau = xout(:,39); test1.lambdar = xout(:,40);

%% test2 无观测器向量场制导
for k=1:Ns
   t = (k-1)*ts;
   % 初始化
   if t==0
       psi = 90*pi/180; x = [0 0 0 10 10 psi]'; tau_c = [0 0]';
       dp = 0.26; B = [1 1;dp -dp];
       Vc = [0 0]'; 
       E = [0 -1;1 0]; I2 = diag([1 1]);
       Tmax = 113.6775; Tmin = -55;
   end
   % 海流变化
%    if t>=50
%        Vc = [0.1 0.6]';
%    end
   if t>=50
       Vc = [V*sin(w*t) V*cos(w*t)]';
   end
   % 外界干扰
   tau_w = [10*sin(0.03*t)*cos(0.2*t)+10 15*cos(0.03*t)*sin(0.2*t)...
            5*cos(0.03*t)*sin(0.2*t)+3]';
   % 无人艇模型
   Tc = B\tau_c;
   [x, xdot, T, fu, fr] = USV2( x, Tc, tau_w, Vc, ts);
   % 计算饱和差值
   tau = B*T;  % 真实值
   Delta_tau = tau_c-tau;  % 饱和差值
   
   % 传感器获得值，及相应计算量
   beta = atan2(x(2),x(1));  % 漂角
   X = x(6)+beta;            % 航向角
   psi = x(6);               % 艏向角
   p = [x(4),x(5)]';             % 无人艇位置向量
   u = x(1); v = x(2); r = x(3); % 无人艇速度和角速度状态
   U = sqrt(x(1)^2+x(2)^2);      % 无人艇总航向速度
   m1 = [cos(psi) sin(psi)]';    % 无人艇姿态向量1
   p_dot = [xdot(4), xdot(5)]';  % 无人艇位置导数向量
   m2 = [-sin(psi) cos(psi)]';   % 无人艇姿态向量2
   vc = v*m2+Vc;                 % 待估计的未知量（即侧滑速度和海流速度）
   % 海流观测器
   [ vchat,vchat_dot,phat] = observer2( p,u,psi,ts );
   
   % 向量场制导
   [rd,e,xd,yd,y,e_alphap,Xd_dot]  = GVFguidance3( pc,ud,p,p_dot,psi,ts );

   % surge controller
   [uhat, fuhat ,tuc, lambdau] = ESOcontroller1( ud,u,tau(1),Delta_tau(1), ts );
   
   % heading controller
   [rhat, frhat, trc, lambdar] = ESOcontroller2( rd,r,tau(2),Delta_tau(2),ts );
   
   % total control signal
   tau_c = [tuc,trc]';

   % 存储时间序列
   xout(k,:) = [x' t e' Vc' vchat' fu fr fuhat frhat tau' tau_c'...
                xd yd ud  rd theta Tc' T' e_alphap Xd_dot vchat_dot' phat' vc'...
                lambdau lambdar];
end

test2.u = xout(:,1); test2.v = xout(:,2); test2.r = xout(:,3);
test2.x = xout(:,4); test2.y = xout(:,5); test2.psi = xout(:,6);
test2.t = xout(:,7); test2.phi1 = xout(:,8); test2.phi2 = xout(:,9);
test2.Vc = xout(:,10:11); test2.vchat = xout(:,12:13);
test2.fu = xout(:,14); test2.fr = xout(:,15);
test2.fuhat = xout(:,16); test2.frhat = xout(:,17);
test2.tau = xout(:,18:19); test2.tau_c = xout(:,20:21);
test2.xd = xout(:,22); test2.yd = xout(:,23);
test2.ud = xout(:,24);
test2.rd = xout(:,25);
test2.theta = xout(:,26);
test2.Tc = xout(:,27:28); test2.T = xout(:,29:30);
test2.e_alphap = xout(:,31); test2.Xd_dot = xout(:,32);
test2.vchat_dot= xout(:,33:34);
test2.xhat = xout(:,35); test2.yhat = xout(:,36);
test2.vc = xout(:,37:38); 
test2.lambdau = xout(:,39); test2.lambdar = xout(:,40);

%% PLOTS

fontsize = 12; fontname = 'TimesNewRoman'; linewid = 2.2; fontweight = 'bold';

figure(1); hold on
xrange=[-40 10 40]; yrange = [-0 10 80];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

plot(test1.yd,test1.xd,'r-','linewid',linewid);
plot(test1.y,test1.x,'g--','linewid',linewid);
plot(test2.y,test2.x,'b--','linewid',linewid);
h1=legend('desired path','trajectory with observer','trajectory without observer');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

% for k=1:4000:Ns
%     pos =[test1.x(k) test1.y(k)]';
%     modelplot(pos,test1.psi(k),xrange,yrange);
% end
% 
% for k=1:4000:Ns
%     pos =[test2.x(k) test2.y(k)]';
%     modelplot1(pos,test2.psi(k),xrange,yrange);
% end
set(gca,'xTick',Xmin:Xinterval:Xmax);
set(gca,'yTick',Ymin:Yinterval:Ymax);
axis([Xmin Xmax,Ymin Ymax]);

h1=xlabel('y (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('x (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

box on
hold off

% 海流估计
figure(2); hold on
plot(test1.t,test1.vc(:,1),'r-',test1.t,test1.vchat(:,1),'b--',...
     test1.t,test1.vc(:,2),'m-',test1.t,test1.vchat(:,2),'g--','linewid',linewid);
h1=legend({'$$n_{3(1)}$$','$$\hat{n}_{3(1)}$$','$$n_{3(2)}$$','$$\hat{n}_{3(2)}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('unknown kinetic disturbances estimation (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(3); hold on
subplot(2,1,1);plot(test1.t,test1.fu,'r-',test1.t,test1.fuhat,'b--','linewid',linewid);
h1=legend('$$f_u$$','$$\hat{f}_u$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$f_{u}\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(test1.t,test1.fr,'r-',test1.t,test1.frhat,'b--','linewid',linewid);
h1=legend('$$f_{r}$$','$$\hat{f}_r$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$f_{r}\ $$(Nm)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(4); hold on
subplot(2,1,1);plot(test1.t,test1.u,'r-',test1.t,test1.ud,'b--','linewid',linewid);
h1=legend('$$U$$','$$U_d$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(test1.t,test1.r,'r-',test1.t,test1.rd,'b--','linewid',linewid);
h1=legend('$$w$$','$$w_d$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(5)
plot(test1.t,test1.theta,'r-','linewid',linewid);
h1=xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\theta$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(6); hold on
plt1 = plot3(test1.x,test1.y,test1.theta,'r-','linewid',linewid);
z = zeros(length(test1.theta),1);
plt2 = plot3(test1.xd,test1.yd, z,'g-','linewid',linewid);
plt3 = plot3(test1.x,test1.y,z,'b--','linewid',linewid);
h1 = legend([plt1,plt2,plt3],{'virtual path','real path','trajetory'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
for k=1:2000:Ns
    plot3([test1.x(k) test1.x(k)],[test1.y(k) test1.y(k)],[test1.theta(k) z(k)],'-m.','linewid',2,'markersize',25);
end
box on

h1=xlabel('x (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('y (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=zlabel('$$\theta\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight); 
box on
hold off

figure(7)
subplot(2,1,1);plot(test1.t,test1.tau(:,1),'linewid',linewid);

h1=ylabel('$$\tau_u\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(test1.t,test1.tau(:,2),'r-','linewid',linewid);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\tau_r\ $$(N.m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(8); hold on
plt1 = plot(test1.t,test1.T(:,1),'r-','linewid',linewid);
plt2 = plot(test1.t,test1.T(:,2),'b-','linewid',linewid);
plt3 = plot(test1.t,test1.Tc(:,1),'g--','linewid',linewid);
plt4 = plot(test1.t,test1.Tc(:,2),'m--','linewid',linewid);
Tmax = repmat(Tmax,length(test1.t),1);
Tmin = repmat(Tmin,length(test1.t),1);
plt5 = plot(test1.t,Tmax,'k--','linewid',linewid);
plt6 = plot(test1.t,Tmin,'c--','linewid',linewid);
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1),plt5(1),plt6(1)],...
         {'$$T_1$$','$$T_2$$','$$T_{c1}$$','$$T_{c2}$$','$$T_{max}$$','$$T_{min}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('Thrust (N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
axes('Position',[0.2,0.6,0.3,0.3]); % 生成子图   
plot(test1.t,test1.T(:,1),'r-',test1.t,test1.T(:,2),'b-',...
     test1.t,test1.Tc(:,1),'g--',test1.t,test1.Tc(:,2),'m--',...
     test1.t,Tmax,'k--',test1.t,Tmin,'c--','linewid',linewid); % 绘制局部曲线图                        
xlim([0,20]); % 设置坐标轴范围  
ylim([-100,150]);
set(gca,'yTick',-100:50:150);
set(gca,'xTick',0:5:20);
box on
hold off

figure(9)
eu = test1.u-test1.ud; er = test1.r-test1.rd;
subplot(3,2,1); plot(test1.t,eu,'r-','linewid',linewid);
h1=ylabel('$$\tilde{u}\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-1 0.2]);
xlim([0 Tfinal]);

subplot(3,2,2); plot(test1.t,er,'r-','linewid',linewid);
h1=ylabel('$$tilde{r}\ $$(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-3 0.5]);
xlim([0 Tfinal]);

subplot(3,2,3); plot(test1.t,test1.e_alphap,'r-','linewid',linewid);
h1=ylabel('$$e_{\alpha p}\ $$(rad)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-2 0.8]);
xlim([0 Tfinal]);

subplot(3,2,4); plot(test1.t,test1.lambdau,'r-','linewid',linewid);
h1=ylabel('$$\lambda_u\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-0.6 0.2]);
xlim([0 Tfinal]);

subplot(3,2,5); plot(test1.t,test1.lambdar,'r-','linewid',linewid);
h1=ylabel('$$\lambda_r\ $$(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-2 0.2]);
xlim([0 Tfinal]);

subplot(3,2,6); plot(test1.t,test1.phi1,'r-',test1.t,test1.phi2,'b-',test1.t,sqrt(test1.phi1.^2+test1.phi2.^2),'g--','linewid',linewid);
h1=legend('$$\phi_1$$','$$\phi_2$$','$$\bar{\phi}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$path following errors $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-0.05 0.15]);
xlim([0 Tfinal]);
box on

figure(10)
efu = test1.fu-test1.fuhat; efr = test1.fr-test1.frhat;
evx = test1.vc(:,1)-test1.vchat(:,1);
evy = test1.vc(:,2)-test1.vchat(:,2);
subplot(221); plot(test1.t,efu,'r-','linewid',linewid);
h1=ylabel('$$e_{fu}\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

subplot(222); plot(test1.t,efr,'r-','linewid',linewid);
h1=ylabel('$$e_{fr}\ $$(Nm)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

subplot(223); plot(test1.t,evx,'r-','linewid',linewid);
h1=ylabel('$$e_{vx}\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

subplot(224); plot(test1.t,evy,'r-','linewid',linewid);
h1=ylabel('$$e_{vy}\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

figure(11)  % 两种制导方法对比

test1.eu = test1.u-test1.ud; test1.er = test1.r-test1.rd;
test2.eu = test2.u-test2.ud; test2.er = test2.r-test2.rd;
subplot(3,2,1); plot(test1.t,test1.eu,'r-',test2.t,test2.eu,'b--','linewid',linewid);
h1=legend('with observer','without observer');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\tilde{u}\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-1 0.2]);
xlim([0 Tfinal]);

subplot(3,2,2); plot(test1.t,test1.er,'r-',test2.t,test2.er,'b--','linewid',linewid);

h1=ylabel('$$\tilde{r}\ $$(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-3 0.5]);
xlim([0 Tfinal]);

subplot(3,2,3); plot(test1.t,test1.e_alphap,'r-',test2.t,test2.e_alphap,'b--','linewid',linewid);

h1=ylabel('$$\tilde{\chi}\ $$(rad)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-2 0.8]);
xlim([0 Tfinal]);

subplot(3,2,4); plot(test1.t,test1.lambdau,'r-',test2.t,test2.lambdau,'b--','linewid',linewid);

h1=ylabel('$$\lambda_u\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-0.6 0.2]);
xlim([0 Tfinal]);

subplot(3,2,5); plot(test1.t,test1.lambdar,'r-',test2.t,test2.lambdar,'b--','linewid',linewid);

h1=ylabel('$$\lambda_r\ $$(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-2 0.2]);
xlim([0 Tfinal]);

subplot(3,2,6); plot(test1.t,sqrt(test1.phi1.^2+test1.phi2.^2),'r-',test2.t,sqrt(test2.phi1.^2+test2.phi2.^2),'b--','linewid',linewid);

h1=ylabel('$$\bar{\phi}\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% ylim([-0.005 0.01]);
xlim([0 Tfinal]);
box on







