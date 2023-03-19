% date: 9 29 2022
% Author: Quyinsong

% 仿真实验1 基于观测器的向量场制导实验
clc
clear all
close all

ts = 0.02;
Tfinal = 250;
Ns = Tfinal/ts;

for k=1:Ns
   t = (k-1)*ts;
   % 初始化
   if t==0
       psi = 90*pi/180; x = [0 0 0 10 10 psi]'; tau_c = [0 0]';
       ud = 0.6; dp = 0.26; B = [1 1;dp -dp];
       Vc = [0.1 0.25]'; 
       E = [0 -1;1 0]; I2 = diag([1 1]);
       Tmax = 113.6775; Tmin = -55;
       pc = 2;
   end
   % 海流变化
   if t>=50
       Vc = [0.15 0.3]';
   end
   if t>=120
       Vc = [0.1*sin(0.05*t) 0.1*cos(0.05*t)]';
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
   [uhat, fuhat ,tuc, lambdau] = NNcontroller1( ud,u,tau(1),Delta_tau(1), ts );
   
   % heading controller
   [rhat, frhat, trc, lambdar] = NNcontroller2( rd,r,tau(2),Delta_tau(2),ts );
   
   % total control signal
   tau_c = [tuc,trc]';

   % 存储时间序列
   xout(k,:) = [x' t e' Vc' vchat' fu fr fuhat frhat tau' tau_c'...
                xd yd ud  rd theta Tc' T' e_alphap Xd_dot vchat_dot' phat' vc'...
                lambdau lambdar];
end

u = xout(:,1); v = xout(:,2); r = xout(:,3);
x = xout(:,4); y = xout(:,5); psi = xout(:,6);
t = xout(:,7); phi1 = xout(:,8); phi2 = xout(:,9);
Vc = xout(:,10:11); vchat = xout(:,12:13);
fu = xout(:,14); fr = xout(:,15);
fuhat = xout(:,16); frhat = xout(:,17);
tau = xout(:,18:19); tau_c = xout(:,20:21);
xd = xout(:,22); yd = xout(:,23);
ud = xout(:,24);
rd = xout(:,25);
theta = xout(:,26);
Tc = xout(:,27:28); T = xout(:,29:30);
e_alphap = xout(:,31); Xd_dot = xout(:,32);
vchat_dot= xout(:,33:34);
xhat = xout(:,35); yhat = xout(:,36);
vc = xout(:,37:38); 
lambdau = xout(:,39); lambdar = xout(:,40);

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
plot(t,vc(:,1),'r-',t,vchat(:,1),'b--',...
     t,vc(:,2),'m-',t,vchat(:,2),'g--','linewid',2);
h1=legend({'$$V_{cx}$$','$$\hat{V}_{cx}$$','$$V_{cy}$$','$$\hat{V}_{cy}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('current estimation (m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(3); hold on
subplot(2,1,1);plot(t,fu,'r-',t,fuhat,'b--','linewid',2);
h1=legend('$$f_U$$','$$\hat{f}_U$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,fr,'r-',t,frhat,'b--','linewid',2);
h1=legend('$$f_{w}$$','$$\hat{f}_w$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(4); hold on
subplot(2,1,1);plot(t,u,'r-',t,ud,'b--','linewid',2);
h1=legend('$$U$$','$$U_d$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,r,'r-',t,rd,'b--','linewid',2);
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
plot(t,theta,'r-','linewid',2);
h1=xlabel('t(s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\theta$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(6); hold on
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

figure(7)
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

figure(8); hold on
plt1 = plot(t,T(:,1),'r-','linewid',2);
plt2 = plot(t,T(:,2),'b-','linewid',2);
plt3 = plot(t,Tc(:,1),'g--','linewid',2);
plt4 = plot(t,Tc(:,2),'m--','linewid',2);
Tmax = repmat(Tmax,length(t),1);
Tmin = repmat(Tmin,length(t),1);
plt5 = plot(t,Tmax,'k--','linewid',2);
plt6 = plot(t,Tmin,'c--','linewid',2);
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1)],{'$$T_1$$','$$T_2$$','$$T_{c2}$$','$$T_{c1}$$','$$T_{max}$$','$$T_{min}$$'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$T\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(9)
eu = u-ud; er = r-rd;
subplot(3,2,1); plot(t,eu,'r-','linewid',2);
h1=ylabel('$$e_u\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-1 0.2]);
xlim([0 Tfinal]);

subplot(3,2,2); plot(t,er,'r-','linewid',2);
h1=ylabel('$$e_r\ $$(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-3 0.5]);
xlim([0 Tfinal]);

subplot(3,2,3); plot(t,e_alphap,'r-','linewid',2);
h1=ylabel('$$e_{\alpha p}\ $$(rad)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-2 0.8]);
xlim([0 Tfinal]);

subplot(3,2,4); plot(t,lambdau,'r-','linewid',2);
h1=ylabel('$$\lambda_u\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-0.6 0.2]);
xlim([0 Tfinal]);

subplot(3,2,5); plot(t,lambdar,'r-','linewid',2);
h1=ylabel('$$\lambda_r\ $$(rad/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-2 0.2]);
xlim([0 Tfinal]);

subplot(3,2,6); plot(t,phi1,'r-',t,phi2,'b-',t,sqrt(phi1.^2+phi2.^2),'g--','linewid',2);
h1=legend('$$\phi_1$$','$$\phi_2$$','$$\bar{\phi}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$\phi\ $$(m)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
ylim([-0.05 0.15]);
xlim([0 Tfinal]);
box on

figure(10)
efu = fu-fuhat; efr = fr-frhat;
evx = vc(:,1)-vchat(:,1);
evy = vc(:,2)-vchat(:,2);
subplot(221); plot(t,efu,'r-','linewid',2);
h1=ylabel('$$e_{fu}\ $$(N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

subplot(222); plot(t,efr,'r-','linewid',2);
h1=ylabel('$$e_{fr}\ $$(Nm)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

subplot(223); plot(t,evx,'r-','linewid',2);
h1=ylabel('$$e_{vx}\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);

subplot(224); plot(t,evy,'r-','linewid',2);
h1=ylabel('$$e_{vy}\ $$(m/s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('time (seconds)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
xlim([0 Tfinal]);










