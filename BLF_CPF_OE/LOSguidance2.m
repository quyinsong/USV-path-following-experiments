function [xd, yd, alpha_psi, psip, kc, up, alpha_r, alpha_u, xe, ye, w, a, b, psie, thetahatout] = LOSguidance2( pc, delta, eta, nu, vs, w, ts, t, ew )
%LOS2 curved path LOS guidace law
%[xd, yd, psid, psif, xe, ye, w] = LOS2( pc, k1, delta, x, y, psi, U, beta, w, ts )
% Author : Quyinsong
% Date: 2022 6 1
% reference : Path-Following Algorithms and Experiments for an Unmanned Surface Vehicle
%%%% INPUT:
% pc : chose curved path
% -------------------pc-----------------------------------
% 1: stright line path; 2: sinusoidal path; 3: circle path
% -------------------pc-----------------------------------
% k1 : k1 make xe converge to 0
% delta : look-ahead distance 
% x,y,psi : USV position and heading angle in north-east reference frame
% U : USV course speed, U = sqrt(u^2+v^2)
% beta : USV sideslip angle, beta = atan2(v,u)
% w : path parameter, xd = xd(w), yd = yd(w)
% ts : sample time , be used to update path parameter w
% t : simulation time, decide when to change path
%%%% OUTPUT:
% xd,yd : path point, xd = xd(w), yd = yd(w)
% psid : disaired course angle
% psif : path tangent-angle to north-axis
% xe,ye : longitude and latitude error described in SF reference frame
% w : path parameter, xd = xd(w), yd = yd(w)

persistent kxe kye kew kr thetahat wx wy
if isempty(kr)
   kxe = 1;
   kye = 1;
   kew = 0.5;
   kr = 5;
   thetahat = 0;
   wx = 0;
   wy = 0;
end
% inputs
x = eta(1); y = eta(2); psi= eta(3);
u = nu(1); v = nu(2); r = nu(3);
U = sqrt(u^2+v^2);
beta = atan2(v,u);

% 过滤beta
persistent betaf;
if isempty(betaf)
   betaf = beta; 
end
betaf_dot = -(betaf-beta)/1.5;
betaf = betaf_dot*ts+betaf;

% 误差约束
vab = 0.05;
sa = 1;
sb = 0.5;
a0 = 15;
b0 = 3;
a = a0*exp(-vab*t)+sa;
b = b0*exp(-vab*t)+sb;
adot = -a0*vab*exp(-vab*t);
bdot = -b0*vab*exp(-vab*t);

% curved path
switch pc
    case 1
        %---------------直线路径-----------------
        xd = 10+5;
        yd = w;
        xd_dw = 0;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
    case 2
        %---------------正弦曲线路径--------------
        k_case2 = 0.5;
        xd = 8*cos(k_case2*w)+8;
        yd = 5*w;
        xd_dw = -8*k_case2*sin(k_case2*w);
        yd_dw = 5;
        xd_ddw = -8*k_case2^2*cos(k_case2*w); yd_ddw = 0;
    case 3
        %---------------圆形路径------------------
        R = 10+5;
        xd = R*cos(w); 
        yd = R*sin(w); 
        xd_dw = -R*sin(w); 
        yd_dw = R*cos(w);
        xd_ddw = -R*cos(w); yd_ddw = -R*sin(w);
    case 4
        %---------------混合路径------------------
        xd = 10+10;
        yd = w;
        xd_dw = 0; xd_ddw = 0;
        yd_dw = 1; yd_ddw = 0;
        if w >= 80 
            R = 10+10;
            xd = R*cos(w-80);
            yd = R*sin(w-80)+80;
            xd_dw = -R*sin(w-80); xd_ddw = -R*cos(w-80);
            yd_dw = R*cos(w-80); yd_ddw = -R*sin(w-80);
        end
        if w >= 80+pi
            xd = -10-10;
            yd = -(w-80-pi)+80;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = -1; yd_ddw = 0;
        end
        if w >= 160+pi
            xd = R*cos(w-160-pi+pi);
            yd = R*sin(w-160-pi+pi);
            xd_dw = -R*sin(w-160-pi+pi); xd_ddw = -R*cos(w-160-pi+pi);
            yd_dw = R*cos(w-160-pi+pi); yd_ddw = -R*sin(w-160-pi+pi);
        end
        if w >= 160+2*pi
            xd = 10+10;
            yd = w-160-2*pi;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = 1; yd_ddw = 0;
        end
end

%---------------------------------------------------------
kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; % 曲线路径的曲率
kcs = (xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3;

psip = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
  
xe = cos(psip) * (x - xd) + sin(psip) * (y - yd);  % USV与虚拟参考点的纵向误差
ye = -sin(psip) * (x - xd) + cos(psip) * (y - yd);  % USV与虚拟参考点的横向误差

% 路径参数更新
wdot =  vs-kew*ew;
w = euler2(wdot, w, ts);

% 制导率设计
ud1 = sqrt(xd_dw^2 + yd_dw^2);
up = wdot*ud1;
ud2 = up*(1-kc*ye*(1-cos(pi*xe^2/2/a^2)^2/cos(pi*ye^2/2/b^2)^2));

Ud = sqrt((ud1*vs)^2+v^2);
kxe0 = sqrt((adot/a)^2+0.5);
rhoe = bdot/b*Ud;
alphae = kye*b^2*sin(pi*(ye+eps)^2/2/b^2)*cos(pi*ye^2/2/b^2)/pi/(ye+eps);
kye0 = (rhoe^2*alphae*ye+rhoe*sqrt(delta*(1-rhoe^2*ye^2)+alphae^2))/(1-rhoe^2*ye^2);
rhoy = kye*b^2*sin(pi*(ye+eps)^2/2/b^2)*cos(pi*ye^2/2/b^2)/pi/(ye+eps)+kye0*ye+thetahat;
rhox = kxe0*xe+kxe*a^2*sin(pi*(xe+eps)^2/2/a^2)*...
          cos(pi*xe^2/2/a^2)/pi/(xe+eps)+thetahat*wx;
alpha_u = 2*u*sin((psi-psip)/2)^2+ud2+v*sin((psi-psip)/2)-rhox;
alpha_psi = psip-betaf-atan2(rhoy,delta);
persistent alpha_psif Tf
if isempty(Tf)
   Tf = 0.1;
   alpha_psif = alpha_psi;
end
alpha_psif_dot = -ssa(alpha_psif-alpha_psi)/Tf;
alpha_psif = ts*alpha_psif_dot+alpha_psif;
psie = ssa(psi-alpha_psif);
qpsi = ssa(alpha_psif-alpha_psi);
rho = (cos(alpha_psi-psip+beta)*sin(psie+qpsi+eps)+sin(alpha_psi-psip+beta)*(cos(psie+qpsi+eps)-1))/(psie+qpsi+eps);
alpha_r = -kr*psie+alpha_psif_dot-(ye/cos(pi*ye^2/2/b^2)^2)*U*rho;

% 自适应律参数更新
Uchat = sqrt(alpha_u^2+v^2);
wx = 0.02*sin(xe)*cos(ye);
wy = Uchat/sqrt(rhoy^2+delta^2);
thetahat_dot = 0.005*(wx*xe/cos(pi*xe^2/2/b^2)^2+wy*ye/cos(pi*ye^2/2/b^2)^2);
thetahat = thetahat_dot*ts+thetahat;
thetahatout = thetahat;

end

