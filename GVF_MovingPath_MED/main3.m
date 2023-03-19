% date: 8 28 2022
% Author: Quyinsong
% reference: 1 Guiding vector field algorithm for a moving
% path following problem
% 2 Guiding vector fields for path following occluded paths
% 同时考虑路径参考坐标系匀速水平移动和动态避障问题
% TimeVaring_nGVF_5 的改动版本，考虑多个障碍物

% 仿真实验1 目标移动 圆形路径 无障碍物
clc
clear all
close all

ts = 0.02;
Tfinal = 500;
Ns = Tfinal/ts;

for k=1:Ns
   t = (k-1)*ts;
   % Initialize
   if t==0
       psi = 30*pi/180; x = [0 0 0 5 -2 psi]'; tau_c = [0 0]';
       Ud = 0.8; dp = 0.26; B = [1 1;dp -dp];
       xT = 0; yT = 0; alphaT = 0*pi/180;
       kr = 0.05; kro = 0.005; vT = 0.3;
       Vc = [0 0]'; 
       RC = 6; E = [0 -1;1 0]; I2 = diag([1 1]);
       Tmax = 113.6775; Tmin = -55;
   end
   tau_w = [10*sin(0.03*t)*cos(0.2*t)+10 15*cos(0.03*t)*sin(0.2*t) 5*cos(0.03*t)*sin(0.2*t)+3]';
   % plant
   Tc = B\tau_c;
   [x, xdot, Tc_1, fU, fw] = USV( x, Tc, tau_w, Vc, ts);
   % 饱和函数近似
   Tc_1a = zeros(2,1);
   for i=1:2
       if Tc(i)<0
           Tc_1a(i) = Tmin*tanh(Tc(i)/Tmin);
       else
           Tc_1a(i) = Tmax*tanh(Tc(i)/Tmax);
       end
   end
   % 计算饱和差值（近似）
   tau_1a = B*Tc_1a;
   Delta_tau_a = tau_c-tau_1a;
   % 计算饱和差值（非近似）
   tau_1 = B*Tc_1;
   Delta_tau = tau_c-tau_1;
   
   beta = atan2(x(2),x(1));
   alpha = x(6)+beta;
   % desired path
   if t>=130
       alphaT = 90*pi/180;
   end
   if t>=260
       alphaT = 180*pi/180;
   end
   if t>=390
       alphaT = 270*pi/180;
   end
   
%    if t>=35
%        alphaT = -90*pi/180;
%    end
   xT_dot = vT*cos(alphaT); xT = euler2(xT_dot,xT,ts); 
   yT_dot = vT*sin(alphaT); yT = euler2(yT_dot,yT,ts);
   R = [cos(alphaT) -sin(alphaT)
        sin(alphaT) cos(alphaT)];
   p = [x(4),x(5)]'; pT = [xT,yT]';
   r = R'*(p-pT); xr = r(1); yr = r(2);
   
   % GVF
   phi = xr^2+yr^2-RC^2;
   n = [2*xr 2*yr]'; tau = [-n(2) n(1)]';
   H = [2 0; 0 2];
% 
%    phi = yr-20*sin(2*pi*xr/50);
%    n = [ -20*2*pi*cos(2*pi*xr/50)/50 1]'; tau = -[-n(2) n(1)]';
%    H = [20*(2*pi)^2*sin(2*pi*xr/50)/2500 0; 0 0];
   
%    phi = yr-8;
%    n = [ 0 1]'; tau = -[-n(2) n(1)]';
%    H = [0 0; 0 0];

%    a = 3; b = 3.5;
%    phi = xr^4+yr^4-2*b^2*(xr^2-yr^2)+a^4-b^4;
%    n = [ 4*xr^3-2*b^2*2*xr, 4*yr^3+2*b^2*2*yr]'; tau = -[-n(2) n(1)]';
%    H = [12*xr^2-2*b^2*2, 0;...
%        0, 12*yr^2-2*b^2*2];
   
   md0 = tau-kr*phi*n;
   md = md0/norm(md0);
   m = [cos(alpha) sin(alpha)]'; mT = [cos(alphaT) sin(alphaT)]';
   mr0 = R'*(Ud*m+Vc-vT*mT);
   mr = mr0/norm(mr0);
   alphad = atan2(md(2),md(1));
   % 计算md_d，alphad_d
   r_d = mr0;  
   md0_d = (E-kr*phi*I2)*H*r_d-kr*n'*r_d*n;
%    md_d = (I2-md0*md0'/norm(md0)^2)*md0_d/norm(md0);
   md_d = -E*(md*md')*E*md0_d/norm(md0);
   alphad_d = -md'*E*md_d;
   % 计算u path following
   e_alpha = mr'*E*md;
   un = alphad_d-1*mr'*E*md;
   up = un*norm(r_d)/(Ud*mr'*R'*m);
   
   % collision avoidance
   if t==0
%        O =[2 1 2 1; 10 8 2 1; 6 10 2 1; 2 6 1.5 0.5; 7 4 2 1;...
%            10 15 2 1; 4 15 2 1; 20 15 3 2]; % 障碍物的横纵坐标、反应域、拒绝域
%        O = [90 20 8 3; 10 85 8 3; 90 92 6 4];
         O = [10 10 0.1 0.01];
%        O = [60 60 10 4; 25 25 15 5];
       NO = length(O(:,1));
   end
   uc = zeros(2,1); uc01 = zeros(2,1); uc0 = zeros(2,1); 
   uc11 = zeros(2,1); uc1 = zeros(2,1); si = zeros(2,1); zi = zeros(2,1);
   for i = 1:NO
       xoi = O(i,1); yoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
       ci = -(Ro1i-Roi)^2;
       psir = (x(4)-xoi)^2+(x(5)-yoi)^2-Roi^2;
       nr = [2*(x(4)-xoi) 2*(x(5)-yoi)]'; taur = [-nr(2) nr(1)]';
       
       mdr0 = taur-kro*psir*nr;
       mdr = mdr0/norm(mdr0);
       
%        mdr01 = -taur-kro*psir*nr;
%        mdr1 = mdr01/norm(mdr01);

       mc0 = Ud*m+Vc;
       mc = mc0/norm(mc0);
       
       mdr0_d = (E-kro*psir*I2)*H*mc0-kro*nr'*mc0*nr;
       mdr_d = -E*(mdr*mdr')*E*mdr0_d/norm(mdr0);
       alphadr_d = -mdr'*E*mdr_d;
       uc01(i) = alphadr_d-0.5*mc'*E*mdr;
       uc0(i) = uc01(i)*norm(mc0)/(Ud*mc'*m);
       
%        mdr01_d = (E-kr*psir*I2)*H*mc0-kro*nr'*mc0*nr;
%        mdr1_d = -E*(mdr1*mdr1')*E*mdr01_d/norm(mdr01);
%        alphadr1_d = -mdr1'*E*mdr1_d;
%        uc11(i) = alphadr1_d-8*mc'*E*mdr1;
%        uc1(i) = uc11(i)*norm(mc0)/(Ud*mc'*m);
%        
%        % 判断左转还是右转
%        if abs(uc0(i)-up)>abs(uc1(i)-up)
%            uc(i) = uc1(i);
%        else
%            uc(i) = uc0(i);
%        end
       
       uc(i) = uc0(i);
%        uc(i) = uc1(i);

       if psir <= ci
           f1 = 0; f2 = 1;
       elseif psir < 0
           l1 = 10; l2 = 5;  % l1 越大则避障优先级越高 l2 越大则路径跟踪优先级越高
           f1 = exp(l1/(ci-psir)); f2 = exp(l2/psir);
%            f1 = tanh(-l1*(ci-psir)); f2 = tanh(-l2*psir);
       else
           f1 = 1; f2 = 0;
       end
       si(i) = f2/(f1+f2); zi(i) = f1/(f1+f2);
   end
   MULzi = 1; SUMuc = 0;
   for i=1:NO
       MULzi = MULzi*zi(i);
       SUMuc = si(i)*uc(i)+SUMuc;
   end
   
   wd = SUMuc+MULzi*up;

   wd = up;

   % surge controller
   U = sqrt(x(1)^2+x(2)^2);
   [Uhat, fUhat, tuc] = controller1( Ud,U,tau_1(1),Delta_tau(1),ts );
   
   % course angular controller
   if U==0
       U = 0.001;
   end
   % 这里调试时间过久，错误原因是将x(3)写为xdot(3)
   w = x(3)+(xdot(2)*x(1)-xdot(1)*x(2))/U^2; % w = r+beta_dot
   [what, fwhat, trc] = controller2( x,wd,w,tau_1(2),Delta_tau(2),ts );
   
   tau_c = [tuc,trc]';

   xout(k,:) = [x' t xT yT phi wd U Uhat fU fUhat w what fw fwhat Tc' Tc_1' Tc_1a' e_alpha...
              alphaT Ud];
end
% x = xout(:,1); y = xout(:,2); xT = xout(:,3); yT = xout(:,4);
u = xout(:,1); v = xout(:,2); r = xout(:,3);
y = xout(:,4); x = xout(:,5); psi = xout(:,6);
t = xout(:,7);
yT = xout(:,8); xT = xout(:,9);
phi = xout(:,10); wd = xout(:,11); U = xout(:,12);
Uhat = xout(:,13); fU  = xout(:,14); fUhat = xout(:,15); 
w = xout(:,16); what = xout(:,17); fw = xout(:,18); fwhat = xout(:,19);
Tc = xout(:,20:21); Tc_1 = xout(:,22:23); Tc_1a = xout(:,24:25);
e_alpha = xout(:,26); alphaT = xout(:,27); Ud = xout(:,28);

%% PLOTS

fontsize = 14; fontname = 'TimesNewRoman'; linewid = 2; fontweight = 'bold';

figure(1); hold on
xrange=[-10 10 50]; yrange = [-10 10 50];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

plot(xT,yT,'b--','linewid',2);
plot(x,y,'r--','linewid',2);

% plot(xT,yT,'b-','linewid',2);
h1=legend('Target path','enclosing path');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% y_tick = 13.5784; x_tick = 13.5521;
% text(x_tick+1.5,y_tick,'\leftarrow Faults occur','color','r','FontSize',12,'linewid',2.5);
% plot(x_tick,y_tick,'*','MarkerSize',15,'color','r')

for k=1:2000:Ns
    pos =[y(k) x(k)]';
    modelplot(pos,psi(k),xrange,yrange);
end

for k=1:2000:Ns
    pos =[yT(k) xT(k)]';
    modelplot2(pos,alphaT(k),xrange,yrange);
    h=rectangle('Position',[xT(k)-RC,yT(k)-RC,2*RC,2*RC],'Curvature',[1,1],'EdgeColor','m','linewid',1);
    set(h,'linestyle','--');
end

plt1 = plot(0,0,'b--','linewid',2);
plt2 = plot(0,0,'r--','linewid',2);
plt3 = plot(0,0,'m--','linewid',1);

h1=legend([plt1(1),plt2(1),plt3(1)],...
    {'target path','trajectory','enclocsing path'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);set(h1,'Location','SouthEast');

% set(gca,'xTick',Xmin:Xinterval:Xmax);
% set(gca,'yTick',Ymin:Yinterval:Ymax);
% axis([Xmin Xmax,Ymin Ymax]);

h1=xlabel('$$y (m)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('$$x (m)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);

% Ma1 = ['A','B','C','D','E','F','G','H','I','J'];
% Ma2 = ['1','2','3','4','5','6','7','8','9','10'];
% l1 = 0.2; l2 = 0.2;

% for k=1:2000:Ns
%     if k==1
%         count = 1;
%     end
%     plot(x(k),y(k),'ks','linewid',2);
%     text(x(k),y(k)+l1,Ma1(count),'FontSize',12,'color','k','FontWeight','bold');
%     plot(xT(k),yT(k),'bo','linewid',2);
%     text(xT(k)+l2,yT(k),Ma2(count),'FontSize',12,'color','b','FontWeight','bold');
%     h=rectangle('Position',[xT(k)-RC,yT(k)-RC,2*RC,2*RC],'Curvature',[1,1],'EdgeColor','m','linewid',1);
%     set(h,'linestyle','--');
%     count=count+1;
% end

% for i=1:NO
% %     xoi = O(i,1); yoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
%     yoi = O(i,1); xoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
%     h=rectangle('Position',[xoi-Roi,yoi-Roi,2*Roi,2*Roi],'Curvature',[1,1],'EdgeColor','k');
%     set(h,'LineStyle','--','linewid',2);
%     h=rectangle('Position',[xoi-Ro1i,yoi-Ro1i,2*Ro1i,2*Ro1i],'Curvature',[1,1],'EdgeColor','k');
%     set(h,'LineStyle','--','linewid',2);
% end

% set(gcf,'unit','centimeters','position',[3 5 15 15]); % 设置图片大小
% set(gca,'position',[0.1 0.1 0.8 0.8]); % 设置坐标轴位置和长度占比
% set(gca,'xtick',-20:20:400,'ytick',-20:20:400); % 设置坐标轴间隔和范围
% xlim([-20 400]); ylim([-20 400]); % 指定坐标轴范围
% axis([-5 15,-5 15]);

% GIF
% 数据点稀疏化
% lx = length(x);
% for i=1:20:lx
%     if i==1
%         k=1;
%     end
%     x1(k,1) = x(i);
%     y1(k,1) = y(i);
%     xT1(k,1) = xT(i);
%     yT1(k,1) = yT(i);
%     k=k+1;
% end
% linecolor = [1 0 0; 0 1 0];
% MovieXY([x1 xT1],[y1 yT1],0.0001,{'o','*'},linecolor);
% Fun_F2gif(F,'TimeVaring_nGVF_5_1.gif',0.01);    % 生成gif图片
box on
hold off

figure(2); hold on
plot(t,phi,'r-','linewid',2);
h1=ylabel('$$\phi(\xi) (m^{2})$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(3); 
plot(t,wd,'r-',t,w,'b-',t,what,'c--','linewid',2);
h1 = legend('$$w_{d}$$','$$w$$','$$\hat{w}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
h1 = xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);

figure(4);
plot(t,Ud,'r-',t,U,'b-',t,Uhat,'c--','linewid',2);
h1=legend('$$U_{d}$$','$$U$$','$$\hat{U}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
box on

figure(5)
plot(t,w-wd,'r-','linewid',2);
h1=ylabel('$$e_{w} (rad/s)$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on

figure(6);hold on
plot(t,Tc(:,1),'r-',t,Tc(:,2),'g-','linewid',2);
plot(t,Tc_1(:,1),'b--',t,Tc_1(:,2),'m--','linewid',2);
h1=legend('$$Tc_{1}$$','$$Tc_{2}$$','$$T_{1}$$','$$T_{2}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
line([0,Tfinal],[Tmax,Tmax],'color',[1 0 0],'linestyle','--','linewid',2);
line([0,Tfinal],[Tmin,Tmin],'color',[1 0 0],'linestyle','--','linewid',2);
h1=xlabel('t (s)'); 
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1=ylabel('Force (N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off


figure(7);hold on
subplot(2,1,1);plot(t,fU,'r-',t,fUhat,'b--','linewid',2);
h1=legend('$$f_{U}$$','$$\hat{f}_{U}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
subplot(2,1,2);plot(t,fw,'r-',t,fwhat,'b--','linewid',2);
h1=legend('$$f_{w}$$','$$\hat{f}_{w}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(8);hold on
plot(t,Tc_1(:,1),'r-',t,Tc_1(:,2),'b-','linewid',2);
h1=legend('$$T_{1}$$','$$T_{2}$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
line([0,Tfinal],[Tmax,Tmax],'color',[1 0 0],'linestyle','--','linewid',2);
line([0,Tfinal],[Tmin,Tmin],'color',[1 0 0],'linestyle','--','linewid',2);
h1 = xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = ylabel('Force (N)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on
hold off

figure(9);
plot(t,e_alpha,'-','linewid',2);
h1 = ylabel('$$sin(e_{\alpha})$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
h1 = xlabel('t (s)');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
box on



