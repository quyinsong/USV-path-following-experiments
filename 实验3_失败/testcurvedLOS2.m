% curved path following test
% date: 6th March 2022
% Author: quyinsong
% Reference: Handbook of marine craft hydrodynamics and motion control
% Edition 2 updated from testcurvedLOS1.m: 4th March 2022 
% 更新部分：使用USV03模型，考虑海流的干扰
% 假设： USV因海流干扰产生的漂角可以精确测得
clc;
clear;
close all;
%% USV parameters
% -------------------惯性参数-----------------------------
m11 = 215; m22 = 265; m33 = 80;
Xu = 70; Xuu = 100; Yv = 100; Yvv = 200; Nr = 50; Nrr = 100;   
%% initial
ts = 0.01;
tfinal = 50;
Ns = tfinal/ts;

psaid_1 = 0;
rd_1 = 0;
xxd_1 = 0;
u_e_sum = 0;
beta_1 = 0;
w = 0;
% USV 
Ustate = [0 0 0 10 10 60*pi/180]';
% kinematic
%% simulation
for k=1:1:Ns
   tout(k,1)=(k-1)*ts;
   % curved path
   %---------------正弦曲线路径-----------------------------
%    K2 = 1;  % 路径虚拟参考点速度更新参数
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
   %---------------直线路径----------------------------------
   ku = 2;  % 路径虚拟参考点速度更新参数
   xd = 30;
   yd = 5*w;
   xd_dw = 0;
   yd_dw = 5;
%    if k*ts >= 20
%       yd_dw = -5; 
%    end
   xd_ddw = 0;
   yd_ddw = 0;
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
   u = Ustate(1); v = Ustate(2); r = Ustate(3);   psai = Ustate(6);
   s = cos(psaif)*(Ustate(4)-xd)+sin(psaif)*(Ustate(5)-yd);  % USV与虚拟参考点的纵向误差
   e = -sin(psaif)*(Ustate(4)-xd)+cos(psaif)*(Ustate(5)-yd);  % USV与虚拟参考点的横向误差
   U = sqrt(u^2+v^2);
   deta = 100;
   beta = atan2(v,u);
   xxsf = psaif-psai-beta;
   Ud = ku*s+U*cos(xxsf);% 路径虚拟参考点的切向速度
   w =  ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  % 路径参数的更新
%    if  k*ts < 20
%       w =  ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  % 路径参数的更新
%    end  
%    if  k*ts >= 20
%       w =  -ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w; % 路径参数的更新，当时间大于给定值时，USV沿着原路径返回
%    end
   xxd = atan2(e,deta);
   xxd_dot = (xxd-xxd_1)/ts; xxd_1 = xxd;   
   beta_dot = (beta-beta_1)/ts; beta_1 = beta;
   kr1 = 1;
   rd = kc*Ud-beta_dot+kr1*(xxsf-xxd)-xxd_dot;
%    rd = psaid_dot-k3*sig(psai_e,m)-k4*sig(psai_e,n); 
   rddot = (rd-rd_1)/ts;
   rd_1 = rd;
   % ———————动力学控制器—————————————   
    %--------------艏向转速控制器------------------------
    r_e = r - rd;
%     Tr =(m22-m11)*u*v+Nr*r+Nrr*abs(r)*r+m33*rddot-m33*k5*sig(r_e,m)-m33*k6*sig(r_e,n);
    kr2 = 0.5;
    fr = (m11-m22)*u*v/m33-Nr*r/m33-Nrr*abs(r)*r/m33;
    Tr = m33*(-fr+rddot-kr2*r_e);% rdot = (m11-m22)*u*v/m33-Nr*r/m33-Nrr*abs(r)*r/m33+Tr/m33+d(3)
    % ———————纵向速度控制器—————————————
    ud = 1;
    u_e = u-ud;
    u_e_sum = u_e+u_e_sum;
    Kpu = 20; Kiu = 0;
    fu = m22*v*r/m11-Xu*u/m11-Xuu*abs(u)*u/m11;
    Fu = m11*(-fu-Kpu*u_e-Kiu*u_e_sum); % udot =  m22*v*r/m11-Xu*u/m11-Xuu*abs(u)*u/m11+Fu/m11+d(1)
%      Fu = 10^5;
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
   Ustate = USV03(Ustate,tao,d);
   % out
   Ustateout(k,:) = Ustate';
   pathout(k,:) = [xd yd];
   seout(k,:) = [s e xxsf-xxd];
   taoout(k,:) = tao';
   nout(k,:) = Ustate(4:6)';
   vout(k,:) = Ustate(1:3)';
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
plot(tout,taoout(:,1),'r-',tout,taoout(:,2),'g-',tout,taoout(:,3),'b-','linewidth',2);
title('控制力和力矩');
xlabel('time/s');ylabel('力/N');
legend('Tx','Ty','Tr');



