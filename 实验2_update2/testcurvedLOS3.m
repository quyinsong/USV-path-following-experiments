% curved path following test
% date: 15st Feb 2022
% Author: quyinsong
% Reference: Handbook of marine craft hydrodynamics and motion control
% Edition 2 updated from testcurvedLOS1.m: 4th March 2022 
% 更新部分：使用USV02模型，考虑海流的干扰
% 假设： USV因海流干扰产生的漂角不可测得
% 处理方法，采用ALOS制导方法，自适应补偿USV产生的漂角
% 参考文献：1.Integral LOS Control for Path Following of Underactuated Marine Surface Vessels in the Presence of Constant Ocean Currents
% 2. Line-of-Sight Path Following for Dubins Paths With Adaptive Sideslip Compensation of Drift Forces
%  实验总结：Fossen 提出的漂角自适应LOS制导率能很好的预测漂角，做出相应补偿，在相对严峻的海流干扰下仍然能保持良好的
%  航向角，跟踪期望路径，实验遇到的问题：USV真实漂角的计算应当如何计算，计算真实漂角是否应当包含海流的速度？
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
% USV 
Ustate = [0 0 0 35 0 pi/2]';
w = 0;
% kinematic
xxd_1 = 0;
beta_hat = 0;
psaif_1 = 0;
rd_1 = 0;
yint = 0;

d = 0*randn(3,1);  % USV动力学干扰
Vc = 1; betac = 0*pi/180; % 海流
current = [Vc betac]';
Vc_I = [Vc*cos(betac) Vc*sin(betac) 0]';

%% simulation
for k=1:1:Ns
   tout(k,1)=(k-1)*ts;
   % curved path
   %---------------正弦曲线路径-----------------------------
   K2 = 1;  % 路径虚拟参考点速度更新参数
   xd = 10*cos(w)+30;
   yd = 5*w;
   xd_dw = -10*sin(w);
   yd_dw = 5;
   xd_ddw = -10*cos(w);
   yd_ddw = 0;
%    if k*ts >= 50
%        xd_dw = 10*sin(w);
%        yd_dw = -5;
%        xd_ddw = 10*cos(w);
%        yd_ddw = 0; 
%    end
   %---------------直线路径----------------------------------
%    K2 = 2;  % 路径虚拟参考点速度更新参数
%    xd = 30;
%    yd = 5*w;
%    xd_dw = 0;
%    yd_dw = 5;
%    xd_ddw = 0;
%    yd_ddw = 0;
%    if k*ts >= 20
%       yd_dw = -5; 
%    end
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
   u = Ustate(1); v = Ustate(2); r = Ustate(3); psai = Ustate(6); 
   s = cos(psaif)*(Ustate(4)-xd)+sin(psaif)*(Ustate(5)-yd);  % USV与虚拟参考点的纵向误差
   e = -sin(psaif)*(Ustate(4)-xd)+cos(psaif)*(Ustate(5)-yd);  % USV与虚拟参考点的横向误差
   % USV的真实漂角----------------------------------
   beta = atan2(v,u);  
   %--------------------------------------------------
   xxsf = psaif-Ustate(6); 
   Ud = sqrt(u^2+v^2)*cos(xxsf)+K2*s;   % 路径虚拟参考点的切向速度
   w =  ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  % 路径参数的更新
%    if  k*ts < 50
%       w =  ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  % 路径参数的更新
%    end  
%    if  k*ts >= 50
%       w =  -ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w; % 路径参数的更新，当时间大于给定值时，USV沿着原路径返回
%    end
   % Kinematic
   deta = 2;
   K1 = 5;
   %--------------LOS艏向导引角----------------------
   gama = 0.03;
   beta_hat_dot = gama*sqrt(u^2+v^2)*deta*e/sqrt(deta^2+(e+deta)^2);
   beta_hat = euler2(beta_hat_dot,beta_hat,ts); % USV的自适应漂角
   xxd = atan2(e+beta_hat*deta,deta);   
   %------------------------------------------
   xxddot = (xxd-xxd_1)/ts;
   xxd_1 = xxd;
   rd = kc*Ud-beta_hat_dot+K1*(xxsf-xxd)-xxddot;  
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
    % -----------------------------------------------------
    fu = (-c13*r-d11*u)/m11;    
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
    %--------------艏向转速控制器------------------------
    Kpr = 5;
    tr =(-Kpr*(r-rd)-fr+rddot)*m0/m22;
    % ———————纵向速度控制器—————————————
    ud = 2;
    eu = u-ud;
    Kpu = 2;
    tu = m11*(-fu-Kpu*eu); % 
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
   tao=[tu 0 tr]';
   Ustate = USV02(Ustate,tao,current,d);
   % out
   Ustateout(k,:) = Ustate';
   pathout(k,:) = [xd yd];
   seout(k,:) = [s e xxsf-xxd];
   taoout(k,:) = tao';
   nout(k,:) = Ustate(4:6)';
   vout(k,:) = Ustate(1:3)';
   betahatout(k) = beta_hat;
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
plot(tout,betaout,'r-',tout,betahatout,'g-','linewidth',2);
title('漂角预测对比');
xlabel('time/s');ylabel('漂角');
legend('实际漂角beta','预测漂角beta_hat');
figure
plot(tout,taoout(:,1),'r-',tout,taoout(:,2),'g-',tout,taoout(:,3),'b-','linewidth',2);
title('控制力和力矩');
xlabel('time/s');ylabel('力/N');
legend('Tx','Ty','Tr');



