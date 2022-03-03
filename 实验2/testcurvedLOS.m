% curved path following test
% date: 15st Feb 2022
% Author: quyinsong
% Reference: Handbook of marine craft hydrodynamics and motion control
% Edition 1 : 12nd Feb 2022 
clc;
clear;
close all;
%% initial
ts = 0.01;
tfinal = 100;
Ns = tfinal/ts;
% USV 
Ustate = [0 0 0 25 5 0]';
w = 0;
% kinematic
xxd_1 = 0;
beta_1 = 0;
psaif_1 = 0;
rd_1 = 0;

%% simulation
for k=1:1:Ns
   tout(k,1)=(k-1)*ts;
   % curved path
   K2 = 30;
   xd = 30*cos(w)+30;
   yd = 30*w;
   xd_dw = -30*sin(w);
   yd_dw = 30;
   xd_ddw = -30*cos(w);
   yd_ddw = 0;
   kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; 
   psaif = atan2(yd_dw,xd_dw);
  
   s = cos(psaif)*(Ustate(4)-xd)+sin(psaif)*(Ustate(5)-yd);
   e = -sin(psaif)*(Ustate(4)-xd)+cos(psaif)*(Ustate(5)-yd);
   beta = atan2(Ustate(2),Ustate(1));
   xxsf = psaif-Ustate(6)-beta;
   Ud = sqrt(Ustate(1)^2+Ustate(2)^2)*cos(xxsf)+K2*s;
   w = ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  
   % Kinematic
   deta = 1;
   K1 = 1.5;
   xxd = atan2(e,deta);
   xxddot = (xxd-xxd_1)/ts;
   xxd_1 = xxd;
   betadot = (beta-beta_1)/ts;
   beta_1 = beta;
   rd = kc*Ud-betadot+K1*(xxsf-xxd)-xxddot;
   rddot = (rd-rd_1)/ts;
   rd_1 = rd;
   % kinetic
    m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
    Nrdot = -1; xg = 0.046; Yrdot = 0;
   % ----------------------------------------------------
    m11 = m-Xudot; 
    m22 = m-Yvdot;
    m23 = m*xg-Yrdot;
    m32 = m*xg-Nvdot;
    m33 = Iz-Nrdot;
   % -----------------------------------------------------
    u = Ustate(1); v = Ustate(2); r = Ustate(3);
   % -----------------------------------------------------
    Mass=[m11  0    0;
          0   m22 m23;
          0   m32 m33];
    m0 = m22*m33-m23*m32;
   %------------------------------------------------------
    Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
    Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
                         Yvv=-36.47287;        Nvv=3.95645;
                        Yrv=-0.805;           Nrv=0.130;
                        Yvr=-0.845;           Nvr=0.080;
                        Yrr=-3.45;            Nrr=-0.75;               
    % ----------------------------------------------------
    c13 = -m*(xg*r+v); 
    c23 = m*u;
    c31 = -c13; c32 = -c23;
    % -----------------------------------------------------
    d11=-Xu-Xuu*abs(u);
    d22=-Yv-Yvv*abs(v)-Yrv*abs(r);
    d23=-Yr-Yvr*abs(v)-Yrr*abs(r);
    d32=-Nv-Nvv*abs(v)-Nrv*abs(r);
    d33=-Nr-Nvr*abs(v)-Nrr*abs(r);
    %-----------------------------------------------------
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
    %-----------------------------------------------------
   Kp = 1;
   tr =(-Kp*(r-rd)-fr*m0/m22+rddot);
   % USV
   d = 0*randn(3,1);
   current = [0 0]';
   tao = [60 0 tr]';
   Ustatedot = USV01(Ustate,tao,current,d);
   Ustate = euler2(Ustatedot,Ustate,ts);
   % out
   Ustateout(k,:) = Ustate';
   pathout(k,:) = [xd yd];
   seout(k,:) = [s e xxsf-xxd];
   taoout(k,:) = tao';
   nout(k,:) = Ustate(4:6)';
   vout(k,:) = Ustate(1:3)';
end
%% plot
figure
plot(pathout(:,2),pathout(:,1),'b-',Ustateout(:,5),Ustateout(:,4),'r-','linewidth',2);
title('ÇúÏßÂ·¾¶¸ú×ÙÐ§¹ûÍ¼');
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
title('¸ú×ÙÎó²î');
xlabel('time/s');ylabel('Îó²î/m');
legend('×ÝÏòÎó²îs','ºáÏòÎó²îe','½Ç¶ÈÎó²îxxsf~');
figure
plot(tout,taoout(:,1),'r-',tout,taoout(:,2),'g-',tout,taoout(:,3),'b-','linewidth',2);
title('¿ØÖÆÁ¦ºÍÁ¦¾Ø');
xlabel('time/s');ylabel('Á¦/N');
legend('Tx','Ty','Tr');



