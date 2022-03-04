function xnext = USV02( x,tao,current,d )
% USV01 xdot = USV( x,tao,taod ) returns the time derivative of 
% the state vector: x = [ u v r x y psi]'  for USV, where
% INPUT: 
% u=x(1): surge velocity (m/s)
% v=x(2): sway velocity (m/s)
% r=x(3): yaw velocity (rad/s)
% x=x(4): position in {n} x-direction (m)
% y=x(5): position in {n} y-direction (m)
% psai=x(6): yaw angle (rad)
% tao=[tu 0 tr]':
% wind=[Vw betaw]'
% current=[Vc betac]'
% OUTPUT: 
% xdot=[udot vdot rdot xdot ydot psaidot Vcdot]':time derivative of state vector
% CS2: mass = 15kg, L = 1.255m, maximum surge force = 2N, maximum yaw force
% = 1.5N, MCLab: L = 40m, B = 6.5m
% 
% Author: Quyinsong
% Data: 2nd March 2022
% Reference: Improved line-of-sight trajectory tracking control of
% under-actuated AUV subjects to ocean currents and input saturation
% USV01版本的更新，考虑海流的影响，建立相对海流的USV运动模型
% 参考文献：1.《Integral LOS Control for Path Following of Underactuated Marine Surface Vessels in the Presence of Constant Ocean Currents》
% 2. MODE LING, IDENTIFICATION, AND ADAPTIVE MANEUVERING OF CYBERSffiP 11: A COMPLETE DESIGN WITH EXPERIMENTS

ts = 0.01;

% check input and state dimentions
if nargin ~=4,error('input number must be 4!');end
if length(x) ~=6,error('state x number must be 6!');end
if length(tao) ~=3,error('ctr input tao number must be 3!');end
if length(current) ~=2,error('current input tao number must be 2!');end
if length(d) ~=3,error('diturbance taod number must be 3!');end

%% USV state:
u=x(1);
v=x(2);
r=x(3);
V = [u v r]';
psai=x(6);
tu = tao(1);
tr = tao(3);
Rbn=[ cos(psai) -sin(psai) 0;
      sin(psai) cos(psai)  0;
      0          0         1];
%% curremt  
Vc= current(1);
betac= current(2);
Vc_I = [Vc*cos(betac) Vc*sin(betac) 0]'; % 惯性坐标系下海流速度
ucx = Vc_I(1); ucy = Vc_I(2);
Vc_B = Rbn'*Vc_I; % 船体坐标系下海流速度
uc = Vc_B(1);  vc = Vc_B(2); 
Vc_B_dot = [r*vc -r*uc 0]'; 
%% relative speed
ur = u-uc;  vr = v-vc;
Vr = [ur vr r]';
%% diturbance
fdu = d(1);
fdv = d(2);
fdr = d(3);
%% USV parameters 
% 以全局变量的形式给出， 在主程序中进行设置
global m; global Xudot; global Nvdot; global Iz; global Yvdot; global Nrdot; global xg; global Yrdot; 
% ------------------------------------------------------
global Xu; global Xuu; global Yv; global Yr; global Yvv; global Yrv; global Yvr; 
global Yrr; global Nv; global Nr; global Nvv; global Nrv; global Nvr; global Nrr;            
% ----------------------------------------------------
global m11; global m22; global m23; global m32; global m33; global m0;
% -----------------------------------------------------
c13 = -m*(xg*r+v); c23 = m*u;
c31 = -c13; c32 = -c23;             % CRB 参数
%------------------------------------------------------
d11=-Xu-Xuu*abs(ur);
d22=-Yv-Yvv*abs(vr)-Yrv*abs(r);
d23=-Yr-Yvr*abs(vr)-Yrr*abs(r);
d32=-Nv-Nvv*abs(vr)-Nrv*abs(r);
d33=-Nr-Nvr*abs(vr)-Nrr*abs(r);    % D 参数
% d11 = 0.9257;
% d22 = 2.8909;  d23 = -0.2601;
% d32 = d23;     d33 = 0.5;
%% matrix expression
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33]; % M = MA+MRB
MA = [-Xudot   0        0;
       0     -Yvdot  -Yrdot;
       0     -Nvdot  -Nrdot];
CRB=[0     0     c13;
     0     0     c23;
    c31   c32     0 ];
CA = [0         0     Yvdot*vr;
      0         0    -Xudot*ur;
    -Yvdot*vr  Xudot*ur   0 ] ;
C = CA+CRB;
D=  [d11  0     0;
     0    d22  d23;
     0    d32  d33];

% time derivatives
% 参考公式： M*vdot+C*v = MA*vcdot+CA*vc-D*vr+tao+d;

Vdot=M\(-C*V+MA*Vc_B_dot+CA*Vc_B-D*Vr+tao)+d ;
Xdot=Rbn*Vr+Vc_I;
xdot=[Vdot ;Xdot];
xnext = euler2(xdot,x,ts);
%% components expression

% detM2=m22*m33-m23*m32;
% m0=detM2;
% 
% fu = (-c13*r-d11*u)/m11;
% fv = (m23*c31*u+m23*c32*v-m33*c23*r-(m33*d22-m23*d32)*v-(m33*d23-m23*d33)*r-m23*tr)/m0;
% fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
% 
% % time derivatives
% 
% xdot = [fu+tu/m11+fdu;
%         fv+fdv;
%         fr+tr*m22/m0+fdr;
%         u*cos(psai)-v*sin(psai)+Vx;
%         u*sin(psai)+v*cos(psai)+Vy;
%         r];

end





