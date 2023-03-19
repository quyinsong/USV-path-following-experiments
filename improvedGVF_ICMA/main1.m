% date: 9 29 2022
% Author: Quyinsong

% 仿真实验1 基于观测器的向量场制导实验
clc
clear all
close all

ts = 0.02;
Tfinal = 180;
Ns = Tfinal/ts;

for k=1:Ns
   t = (k-1)*ts;
   % 初始化
   if t==0
       psi = 90*pi/180; x = [0 0 0 10 10 psi]'; tau_c = [0 0]';
       Ud = 0.8; dp = 0.26; B = [1 1;dp -dp];
       Vc = [0.1 0.2]'; 
       E = [0 -1;1 0]; I2 = diag([1 1]);
       Tmax = 113.6775; Tmin = -55;
       pc = 2;
   end
   % 海流变化
   if t>=50
       Vc = [0.15 0.25]';
   end
   
   % 外界干扰
   tau_w = [10*sin(0.03*t)*cos(0.2*t)+10 15*cos(0.03*t)*sin(0.2*t)...
            5*cos(0.03*t)*sin(0.2*t)+3]';
   % 无人艇模型
   Tc = B\tau_c;
   [x, xdot, T, fU, fw] = USV( x, Tc, tau_w, Vc, ts);
   % 计算饱和差值
   tau = B*T;  % 真实值
   Delta_tau = tau_c-tau;  % 饱和差值
   
   % 传感器获得值，及相应计算量
   beta = atan2(x(2),x(1));  % 漂角
   X = x(6)+beta;            % 航向角
   p = [x(4),x(5)]';         % 无人艇位置向量
   U = sqrt(x(1)^2+x(2)^2);  % 无人艇总航向速度
   m = [cos(X) sin(X)]';     % 无人艇姿态向量
   p_dot = [xdot(4), xdot(5)]';  %  无人艇位置导数向量
   % 海流观测器
   [ Vchat,Vchat_dot,phat] = observer( p,U,X,ts );
   
   % 向量场制导
   [wd,e,xd,yd,theta,e_alphap,Xd_dot]  = GVFguidance( pc,Ud,Vchat,Vchat_dot,p,p_dot,X,ts );

   % surge controller
   [Uhat, fUhat, tuc] = controller1(Ud,U,tau(1),Delta_tau(1),ts );
   
   % heading controller
   w = x(3)+(xdot(2)*x(1)-xdot(1)*x(2))/U^2; % w = r+beta_dot
   [what, fwhat, trc] = controller2( x,wd,w,tau(2),Delta_tau(2),ts );
   
   % total control signal
   tau_c = [tuc,trc]';

   % 存储时间序列
   xout(k,:) = [x' t e' Vc' Vchat' fU fw fUhat fwhat tau' tau_c'...
                xd yd U Ud w wd theta Tc' T' e_alphap Xd_dot Vchat_dot' phat'];
end

u = xout(:,1); v = xout(:,2); r = xout(:,3);
x = xout(:,4); y = xout(:,5); psi = xout(:,6);
t = xout(:,7); phi1 = xout(:,8); phi2 = xout(:,9);
Vc = xout(:,10:11); Vchat = xout(:,12:13);
fU = xout(:,14); fw = xout(:,15);
fUhat = xout(:,16); fwhat = xout(:,17);
tau = xout(:,18:19); tau_c = xout(:,20:21);
xd = xout(:,22); yd = xout(:,23);
U = xout(:,24); Ud = xout(:,25);
w = xout(:,26); wd = xout(:,27);
theta = xout(:,28);
Tc = xout(:,29:30); T = xout(:,31:32);
e_alphap = xout(:,33); Xd_dot = xout(:,34);
Vchat_dot= xout(:,35:36);
xhat = xout(:,37); yhat = xout(:,38);

%% PLOTS

fontsize = 12; fontname = 'TimesNewRoman'; linewid = 2; fontweight = 'bold';

figure(1); hold on
xrange=[-40 10 40]; yrange = [-0 10 80];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

plot(yd,xd,'b-','linewid',2);
plot(y,x,'m--','linewid',2);
h1=legend('desired path','trajectory');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

for k=1:1000:Ns
    pos =[x(k) y(k)]';
    modelplot(pos,psi(k),xrange,yrange);
end

h1=xlabel('y (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('x (m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

box on
hold off

figure(2); hold on
subplot(2,1,1); plot(t,phi1,'r-',t,phi2,'b-','linewid',2);
h1=legend('$$\phi_1$$','$$\phi_2$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\phi\ (m^2)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2); plot(t,sqrt(phi1.^2+phi2.^2),'r-','linewid',2);
h1=text(90,0.1,'$\bar{\phi}=\sqrt{\phi^{2}_{1}+\phi^{2}_{2}}$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\bar{\phi}\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(3); hold on
plot(t,Vc(:,1),'r-',t,Vchat(:,1),'b--',...
            t,Vc(:,2),'m-',t,Vchat(:,2),'g--','linewid',2);
h1=legend({'$$V_{cx}$$','$$\hat{V}_{cx}$$','$$V_{cy}$$','$$\hat{V}_{cy}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(4); hold on
subplot(2,1,1);plot(t,fU,'r-',t,fUhat,'b--','linewid',2);
h1=legend('$$f_U$$','$$\hat{f}_U$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,fw,'r-',t,fwhat,'b--','linewid',2);
h1=legend('$$f_{w}$$','$$\hat{f}_w$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(5); hold on
subplot(2,1,1);plot(t,U,'r-',t,Ud,'b--','linewid',2);
h1=legend('$$U$$','$$U_d$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,w,'r-',t,wd,'b--','linewid',2);
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

figure(6)
plot(t,theta,'r-','linewid',2);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\theta$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(7); hold on
plt1 = plot3(x,y,theta,'r-','linewid',2);
z = zeros(length(theta),1);
plt2 = plot3(xd,yd,z,'g-','linewid',2);
plt3 = plot3(x,y,z,'b--','linewid',2);
h1 = legend([plt1,plt2,plt3],{'virtual path','real path','trajetory'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
for k=1:2000:Ns
    plot3([x(k) x(k)],[y(k) y(k)],[theta(k) z(k)],'-m.','linewid',2,'markersize',25);
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

figure(8)
subplot(2,1,1);plot(t,tau(:,1),'linewid',2);

h1=ylabel('$$\tau_u\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,tau(:,2),'r-','linewid',2);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\tau_r\ $$(N.m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(9); hold on
plt1 = plot(t,T(:,1),'r-','linewid',2);
plt2 = plot(t,T(:,2),'b-','linewid',2);
Tmax = repmat(Tmax,length(t),1);
Tmin = repmat(Tmin,length(t),1);
plt3 = plot(t,Tmax,'g--','linewid',2);
plt4 = plot(t,Tmin,'m--','linewid',2);
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1)],{'$$T_1$$','$$T_2$$','$$T_{max}$$','$$T_{min}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$T\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(10)
subplot(2,1,1);plot(t,e_alphap,'linewid',2);
h1=ylabel('$$e_{\alpha}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,w-wd,'r-','linewid',2);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$e_{w}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(11)
plot(t,Xd_dot,'linewid',2);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\dot{Xd}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(12)
plot(t,Vchat_dot(:,1),'r-',t,Vchat_dot(:,2),'b--','linewid',2);
h1=legend('$$\dot{\hat{V}}_{cx}$$','$$\dot{\hat{V}}_{cy}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(13)
subplot(2,1,1); plot(t,x,'r-',t,xhat,'b--','linewid',2);
h1=legend('$$x$$','$$\hat{x}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2); plot(t,y,'r-',t,yhat,'b--','linewid',2);
h1=legend('$$y$$','$$\hat{y}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

