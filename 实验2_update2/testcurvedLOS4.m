% curved path following test
% date: 6th March 2022
% Author: quyinsong
% Reference: 1.Handbook of marine craft hydrodynamics and motion control
% 2. 固定时间预测器下的欠驱动无人艇路径跟踪控制
% Edition 2 updated from testcurvedLOS1.m: 4th March 2022 
% 更新部分：使用USV02模型，考虑海流的干扰，使用参考文献2中的漂角观测器
% 假设： USV因海流干扰产生的漂角不能精确测得
% 实验总结： 在海流干扰较小的情况下能很好的预测漂角，但当海流变化时，参数调节困难，路径跟踪效果不好
clc;
clear;
close all;
%% USV parameters
% -------------------惯性参数-----------------------------
global m; global Xudot; global Nvdot; global Iz; global Yvdot; global Nrdot; global xg; global Yrdot;
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
global m11; global m22; global m23; global m32; global m33; global m0;             
m11 = 25.8; m22 = 33.8; m23 = 1.0948; m32 = m23;  m33 = 2.76; m0 = 92.0894; % M 参数
% m11 = m-Xudot; m22 = m-Yvdot; m23 = m*xg-Yrdot;
% m32 = m*xg-Nvdot; m33 = Iz-Nrdot; m0 = m22*m33-m23*m32;  
% -------------------水动力参数---------------------------
global Xu; global Xuu; global Yv; global Yr; global Yvv; global Yrv; global Yvr; 
global Yrr; global Nv; global Nr; global Nvv; global Nrv; global Nvr; global Nrr; 
Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
                     Yvv=-36.47287;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;    
%% initial
ts = 0.01;
tfinal = 100;
Ns = tfinal/ts;

s_hat = 5;
e_hat = 5;
beta_hat = 0.5;

xxd_1 = 0;
rd_1 = 0.5;
beta_1 =0;
w = 0;
% USV 
x = [0 0 0 35 10 0]';
% kinematic
kx1 = 0.5; ky1 = 0.5; kx2 = 0.5; ky2 = 0.5;
k3 = 0.5; k4 = 0.6; k5 = 0.5; k6 = 0.6; 
ku = 2;
mm = 0.3; nn = 1.8;
C = 0.8;
%% simulation
for k=1:1:Ns
   tout(k,1)=(k-1)*ts;
   % curved path
   %---------------正弦曲线路径-----------------------------
%    xd = 10*cos(w)+30;
%    yd = 5*w;
%    xd_dw = -10*sin(w);
%    yd_dw = 5;
%    xd_ddw = -10*cos(w);
%    yd_ddw = 0;
%    if k*ts >= 50
%        xd_dw = 10*sin(w);
%        yd_dw = -5;
%        xd_ddw = 10*cos(w);
%        yd_ddw = 0; 
%    end
   %---------------拟正弦曲线路径-----------------------------
   % 参考文献：固定时间预测器下的欠驱动无人艇路径跟踪控制
%    xd = 10*cos(0.3*w)+3*w+15;
%    yd = 3*w;
%    xd_dw = -3*sin(0.3*w)+3;
%    yd_dw = 3;
%    xd_ddw = -0.9*cos(0.3*w);
%    yd_ddw = 0;
   %---------------直线路径----------------------------------
   xd = 40;
   yd = 5*w;
   xd_dw = 0;
   yd_dw = 5;
   xd_ddw = 0;
   yd_ddw = 0;
% %    if k*ts >= 20
% %       yd_dw = -5; 
% %    end
   %---------------圆形路径----------------------------------
%    K2 = 20;  % 路径虚拟参考点速度更新参数 
%    R = 10;
%    xd = R*cos(w)+20; 
%    yd = R*sin(w)+20; 
%    xd_dw = -10*sin(w); 
%    yd_dw = 10*cos(w);
%    xd_ddw = -10*cos(w); yd_ddw = -10*sin(w);
   % --------------------------------------------------------
   kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; % 曲线路径的曲率
%    kc = 1/R; % 当曲线为圆时的曲率；
   psaif = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
%    psaif = w+pi/2;  % 当曲线为圆时的切向角
   %---------------------运动学控制器------------------------------------    
   u = x(1); v = x(2); r = x(3);   
   s = cos(psaif)*(x(4)-xd)+sin(psaif)*(x(5)-yd);  % USV与虚拟参考点的纵向误差
   e = -sin(psaif)*(x(4)-xd)+cos(psaif)*(x(5)-yd);  % USV与虚拟参考点的横向误差
   U = sqrt(u^2+v^2);
   beta = atan2(v,u);
%    xxsf = psaif-x(6)-beta; 
   xxsf = psaif-x(6)-beta_hat;
   Ud = ku*s_hat+U*cos(xxsf);% 路径虚拟参考点的切向速度
   s_e = s_hat-s;
   e_e = e_hat-e;     
   s_hat_dot = U*cos(x(6)-psaif)-U*sin(x(6)-psaif)*beta_hat+kc*Ud*e_hat-Ud-kx1*sig(s_e,mm)-kx2*sig(s_e,nn);
   s_hat = euler2(s_hat_dot,s_hat,ts);
   e_hat_dot = U*sin(x(6)-psaif)+U*cos(x(6)-psaif)*beta_hat-kc*Ud*s_hat-ky1*sig(e_e,mm)-ky2*sig(e_e,nn);
   e_hat = euler2(e_hat_dot,e_hat,ts);
   beta_hat_dot = C*(U*sin(x(6)-psaif)*s_e-U*cos(x(6)-psaif)*e_e);
   beta_hat = euler2(beta_hat_dot,beta_hat,ts);
   w =  ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  % 路径参数的更新
%    if  k*ts < 20
%       w =  ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  % 路径参数的更新
%    end  
%    if  k*ts >= 20
%       w =  -ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w; % 路径参数的更新，当时间大于给定值时，USV沿着原路径返回
%    end
   % Kinematic
   deta = 2;
   %--------------LOS艏向导引角----------------------
   xxd = atan2(e,deta);
   %------------------------------------------
   xxd_dot = (xxd-xxd_1)/ts;  xxd_1 = xxd;
%    beta_dot = (beta-beta_1)/ts; beta_1 = beta;
   xxd_e = xxsf-xxd;
   k1 = 0;
   rd = kc*Ud-beta_hat_dot-xxd_dot+k3*sig(xxd_e,mm)+k4*sig(xxd_e,nn)+k1*xxd_e; % 当采用传统LOS制导率时，假设漂角beta可以精确测得
   rddot = (rd-rd_1)/ts;
   rd_1 = rd;
   % ———————动力学控制器—————————————
    c13=-m23*r-m22*v; c23 = m11*u;
    c31 = -c13; c32 = -c23;
    d11=-Xu-Xuu*abs(u);
    d22=-Yv-Yvv*abs(v)-Yrv*abs(r);
    d23=-Yr-Yvr*abs(v)-Yrr*abs(r);
    d32=-Nv-Nvv*abs(v)-Nrv*abs(r);
    d33=-Nr-Nvr*abs(v)-Nrr*abs(r);    
    %--------------艏向转速控制器------------------------
    r_e = r - rd;
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
    Kpr = 0;
    Tr =(rddot-k5*sig(r_e,mm)-k6*sig(r_e,nn)-fr-Kpr*r_e)*m0/m22;
    % 论文2中提出的固定时间观测器本质上是增加阶数可自定的阻尼项，当mm=1时，为线性阻尼,当mm=2时，为非线性阻尼项
    % ———————纵向速度控制器—————————————
    ud = 2;
    u_e = u-ud;
    Kpu = 2; 
    fu = (-c13*r-d11*u)/m11;
    Fu = (-Kpu*u_e-fu)*m11; 
%    % ———————纵向推力输出限幅———————————
%    tu_limit = 40;
%    if abs(tu) >= tu_limit
%       tu = sign(tu)*tu_limit;
%    end
%    %---------------转艏力矩输出限幅--------------------------
%    tr_limit = 20;
%    if abs(tr) >= tr_limit
%       tr = sign(tr)*tr_limit;
%    end
   %----------------USV状态更新------------------------------
   tao=[Fu 0 Tr]';
   d = 0*randn(3,1);
   x = USV02(x,tao,[0.5 0]',d);
   % out
   Ustateout(k,:) = x';
   pathout(k,:) = [xd yd];
   seout(k,:) = [s e xxsf-xxd];
   sehatout(k,:) = [s_hat e_hat];
   taoout(k,:) = tao';
   nout(k,:) = x(4:6)';
   vout(k,:) = x(1:3)';
   beta_hatout(k) = beta_hat;
   betaout(k) = beta;
   
end
%% plot
for k=1:1:Ns
    pos =[Ustateout(k,4) Ustateout(k,5)]';
    psai(k)=Ustateout(k,6);
    if k==1
        modelplot(pos,psai(k));
    end
    if rem(k,500)==0
        modelplot(pos,psai(k));
    end   
end
plot(pathout(:,2),pathout(:,1),'b-','linewidth',1);
plot(Ustateout(:,5),Ustateout(:,4),'r--','linewidth',1)
hold off
figure
plot(pathout(:,2),pathout(:,1),'b-',Ustateout(:,5),Ustateout(:,4),'r-','linewidth',2);
title('曲线路径跟踪效果图');
xlabel('E/m');ylabel('N/m');
figure
subplot(3,1,1); plot(tout,Ustateout(:,1),'r-','linewidth',2);title('u');xlabel('time/s');ylabel('u(m/s)');
subplot(3,1,2); plot(tout,Ustateout(:,2),'r-','linewidth',2);title('v');xlabel('time/s');ylabel('v(m/s)');
subplot(3,1,3); plot(tout,Ustateout(:,3),'r-','linewidth',2);title('r');xlabel('time/s');ylabel('r(rad/s)');
figure
subplot(3,1,1); plot(tout,Ustateout(:,4),'r-','linewidth',2);title('x');xlabel('time/s');ylabel('x(m)');
subplot(3,1,2); plot(tout,Ustateout(:,5),'r-','linewidth',2);title('y');xlabel('time/s');ylabel('y(m)');
subplot(3,1,3); plot(tout,Ustateout(:,6)*180/pi,'r-','linewidth',2);title('psai');xlabel('time/s');ylabel('psai(rad)');
figure
plot(tout,seout(:,1),'r-',tout,seout(:,2),'b-',tout,seout(:,3),'k-','linewidth',2);
title('跟踪误差');
xlabel('time/s');ylabel('误差/m');
legend('纵向误差s','横向误差e','角度误差xxsf~');
figure
plot(tout,seout(:,1),'r-',tout,seout(:,2),'b-',tout,sehatout(:,1),'k-',tout,sehatout(:,2),'g-','linewidth',2);
title('跟踪误差');
xlabel('time/s');ylabel('误差/m');
legend('纵向误差s','横向误差e','s_hat','e_hat');
figure
plot(tout,betaout,'r-',tout,beta_hatout,'b-','linewidth',2);
title('侧滑角预测效果');
xlabel('time/s');ylabel('预测值');
legend('实际侧滑角beta','预测侧滑角beta_hat');
figure
plot(tout,taoout(:,1),'r-',tout,taoout(:,2),'g-',tout,taoout(:,3),'b-','linewidth',2);
title('控制力和力矩');
xlabel('time/s');ylabel('力/N');
legend('Tx','Ty','Tr');




