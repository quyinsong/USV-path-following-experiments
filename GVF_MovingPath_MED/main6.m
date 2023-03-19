% date: 8 28 2022
% Author: Quyinsong
% reference: 1 Guiding vector field algorithm for a moving
% path following problem
% 2 Guiding vector fields for path following occluded paths
% 同时考虑路径参考坐标系匀速水平移动和动态避障问题
% TimeVaring_nGVF_5 的改动版本，考虑多个障碍物 

% 仿真实验1 目标移动 圆形路径 有障碍物 且障碍物移动
clc
clear all
close all

ts = 0.02;
Tfinal = 350;
Ns = Tfinal/ts;

for k=1:Ns
   t = (k-1)*ts;
   % Initialize
   if t==0
       psi = 30*pi/180; x = [0 0 0 5 -2 psi]'; tau_c = [0 0]';
       Ud = 0.8; dp = 0.26; B = [1 1;dp -dp];
       xT = 0; yT = 0; alphaT = 45*pi/180; vT = 0.3;
       xO = 35; yO = 30; alphaO = 225*pi/180; vO = 0.1;
       xO1 = 70; yO1 = 30; alphaO1 = 135*pi/180; vO1 = 0.1;
       kr = 0.005; kro = 0.005; 
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
%    if t>=130
%        alphaT = 90*pi/180;
%    end
%    if t>=260
%        alphaT = 180*pi/180;
%    end
%    if t>=390
%        alphaT = 270*pi/180;
%    end
   
%    if t>=35
%        alphaT = -90*pi/180;
%    end
   xT_dot = vT*cos(alphaT); xT = euler2(xT_dot,xT,ts); 
   yT_dot = vT*sin(alphaT); yT = euler2(yT_dot,yT,ts);
   R = [cos(alphaT) -sin(alphaT)
        sin(alphaT) cos(alphaT)];
   p = [x(4),x(5)]'; pT = [xT,yT]';
   r = R'*(p-pT); xr = r(1); yr = r(2);
   
   xO_dot = vO*cos(alphaO); xO = euler2(xO_dot,xO,ts); 
   yO_dot = vO*sin(alphaO); yO = euler2(yO_dot,yO,ts);
   xO1_dot = vO1*cos(alphaO1); xO1 = euler2(xO1_dot,xO1,ts); 
   yO1_dot = vO1*sin(alphaO1); yO1 = euler2(yO1_dot,yO1,ts); 
   % GVF
   phi = xr^2+yr^2-RC^2;
   n = [2*xr 2*yr]'; tau = [-n(2) n(1)]';
   H = [2 0; 0 2];
   
   md0 = tau-kr*phi*n;
   md = md0/norm(md0);
   m = [cos(alpha) sin(alpha)]'; mT = [cos(alphaT) sin(alphaT)]';
   mr0 = R'*(Ud*m+Vc-vT*mT);
   mr = mr0/norm(mr0);
   alphad = atan2(md(2),md(1));
   % 计算md_d，alphad_d
   r_d = mr0;  
   md0_d = (E-kr*phi*I2)*H*r_d-kr*n'*r_d*n;
   md_d = -E*(md*md')*E*md0_d/norm(md0);
   alphad_d = -md'*E*md_d;
   % 计算u path following
   e_alpha = mr'*E*md;
   un = alphad_d-2*mr'*E*md;
   up = un*norm(r_d)/(Ud*mr'*R'*m);
   
   % collision avoidance
   O = [xO yO 8 3 vO alphaO; xO1 yO1 8 3 vO1 alphaO1];
   NO = length(O(:,1));
   uc = zeros(2,1); uco0 = zeros(2,1); uco = zeros(2,1);
   si = zeros(2,1); zi = zeros(2,1);
   for i = 1:NO
       xoi = O(i,1); yoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
       voi = O(i,5); alphaOi = O(i,6);
       ci = -(Ro1i-Roi)^2;
       pOi = [xoi,yoi]';
       Ri = [cos(alphaOi) -sin(alphaOi)
             sin(alphaOi) cos(alphaOi)];
       rO = Ri'*(p-pOi); xro = rO(1); yro = rO(2); 
       psir = xro^2+yro^2-Roi^2;
       no = [2*xro 2*yro]'; tauo = [-no(2) no(1)]';
       
       mdo0 = tauo-kro*psir*no;
       mdo = mdo0/norm(mdo0);
       
       m = [cos(alpha) sin(alpha)]'; mo = [cos(alphaOi) sin(alphaOi)]';
       mro0 = Ri'*(Ud*m+Vc-voi*mo);
       mro = mro0/norm(mro0);
       % cal alphado_d
       alphado = atan2(mro(2),mro(1));
       ro_d = mro0;
       mdo0_d = (E-kro*psir*I2)*H*ro_d-kro*no'*ro_d*no;
       mdo_d = -E*(mdo*mdo')*E*mdo0_d/norm(mdo0);
       alphado_d = -mdo'*E*mdo_d;
       
       e_alphao = mro'*E*mdo;
       uco0(i) = alphado_d-0.5*mro'*E*mdo;
       uco(i) = uco0(i)*norm(ro_d)/(Ud*mro'*Ri'*m);
       
       uc(i) = uco(i);

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
              alphaT Ud xO yO xO1 yO1];
end
% x = xout(:,1); y = xout(:,2); xT = xout(:,3); yT = xout(:,4);
u = xout(:,1); v = xout(:,2); r = xout(:,3);
y = xout(:,4); x = xout(:,5); psi = xout(:,6);
t = xout(:,7);
yT = xout(:,8); xT = xout(:,9);alphaT = xout(:,27);
phi = xout(:,10);
wd = xout(:,11);w = xout(:,16);what = xout(:,17);
yO = xout(:,29); xO = xout(:,30); yO1 = xout(:,31); xO1 = xout(:,32);

%% PLOTS
tfinal1 = 350;
Ns1 = tfinal1/ts; % 截取特定时间段的数据

y = y(1:Ns1); x = x(1:Ns1); psi = psi(1:Ns1);
yT = yT(1:Ns1); xT = xT(1:Ns1);alphaT = alphaT(1:Ns1);
wd = y(1:Ns1);w = y(1:Ns1); what = what(1:Ns1);
yO = yO(1:Ns1); xO = xO(1:Ns1); yO1 = yO1(1:Ns1); xO1 = xO1(1:Ns1);

fontsize = 14; fontname = 'TimesNewRoman'; linewid = 2; fontweight = 'bold';

figure(1); hold on
xrange=[-10 10 100]; yrange = [-10 10 100];

Xmax = xrange(3);  Xinterval = xrange(2); Xmin = xrange(1);
Ymax = yrange(3);  Yinterval = yrange(2); Ymin = yrange(1);

% 时间标注
h1=text(40,90,'$$t=350s$$');
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);
% 绘制障碍物，围捕目标以及本体轨迹
plt1=plot(xO,yO,'g-','linewid',2);

plt2=plot(xO1,yO1,'c-','linewid',2);

plt3=plot(xT,yT,'b-','linewid',2);

plt4=plot(x,y,'r-','linewid',2);

% plot(xT,yT,'b-','linewid',2);

% y_tick = 13.5784; x_tick = 13.5521;
% text(x_tick+1.5,y_tick,'\leftarrow Faults occur','color','r','FontSize',12,'linewid',2.5);
% plot(x_tick,y_tick,'*','MarkerSize',15,'color','r')

% 绘制箭头
ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',10,'HeadWidth',5,'color','g','linewidth',2);
set(ah,'parent',gca);
set(ah,'position',[30 35 -8 -8]);

ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',10,'HeadWidth',5,'color','c','linewidth',2);
set(ah,'parent',gca);
set(ah,'position',[30 70 8 -8]);

ah = annotation('arrow',...
            'headStyle','cback1','HeadLength',10,'HeadWidth',5,'color','r','linewidth',2);
set(ah,'parent',gca);
set(ah,'position',[-2 5 10*sin(30*pi/180) 10*cos(30*pi/180)]);
% 绘制无人艇和待围捕船
for k=1:1000:Ns1
    pos =[y(k) x(k)]';
    modelplot(pos,psi(k),xrange,yrange);
end

for k=1:1000:Ns1
    pos =[yT(k) xT(k)]';
    modelplot2(pos,alphaT(k),xrange,yrange);
    plt5=rectangle('Position',[xT(k)-RC,yT(k)-RC,2*RC,2*RC],'Curvature',[1,1],'EdgeColor','m','linewid',1);
    set(plt5,'linestyle','--');
end

plt5 = plot(0,0,'m--','linewid',1);
h1=legend([plt1(1),plt2(1),plt3(1),plt4(1),plt5(1)],...
    {'Obstacle1','Obstacle2','target path','trajectory','enclocsing path'});
set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
set(h1,'FontWeight',fontweight);set(h1,'Location','SouthEast');

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

O = [xO,xO1,yO,yO1]; coloro = ['g','c'];
for i=1:NO
%     xoi = O(i,1); yoi = O(i,2); Roi = O(i,3); Ro1i = O(i,4);
    yoi = O(Ns1,i+2); xoi = O(Ns1,i); Roi = 8; Ro1i = 3;
    h=rectangle('Position',[xoi-Roi,yoi-Roi,2*Roi,2*Roi],'Curvature',[1,1],'EdgeColor',coloro(i));
    set(h,'LineStyle','-','linewid',2);
    h=rectangle('Position',[xoi-Ro1i,yoi-Ro1i,2*Ro1i,2*Ro1i],'Curvature',[1,1],'EdgeColor',coloro(i));
    set(h,'LineStyle','-','linewid',2);
    
end

% GIF
% 数据点稀疏化
lx = length(x);
for i=1:30:lx
    if i==1
        k=1;
    end
    x1(k,1) = x(i);
    y1(k,1) = y(i);
    xT1(k,1) = xT(i);
    yT1(k,1) = yT(i);
    xOf(k,1) = xO(i);
    yOf(k,1) = yO(i);
    xO1f(k,1) = xO1(i);
    yO1f(k,1) = yO1(i);
    k=k+1;
end
linecolor = [1 0 0; 0 1 0; 0 0 1; 0 1 1];
F=MovieXY([x1 xT1 xOf xO1f],[y1 yT1 yOf yO1f],0.0001,{'o','*','s','+'},linecolor);
Fun_F2gif(F,'TimeVaring_nGVF_5_1.gif',0.01);    % 生成gif图片
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

% figure(3); 
% plot(t,wd,'r-',t,w,'b-',t,what,'c--','linewid',2);
% h1 = legend('$$w_{d}$$','$$w$$','$$\hat{w}$$');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);
% h1 = xlabel('t (s)');
% set(h1,'Interpreter','latex'); set(h1,'FontSize',fontsize); set(h1,'FontName',fontname);








