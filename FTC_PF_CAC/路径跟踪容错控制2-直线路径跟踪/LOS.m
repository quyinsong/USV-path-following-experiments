function [xd, yd, Xsf_d, psif, kc, up , xe, ye, w] = LOS( pc, k1, delta, x, y, psi, U, beta, w, ts )
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

% check input and state dimentions
if nargin ~=10,error('input number must be 10!');end

% curved path
switch pc
    case 1
        %---------------直线路径-----------------
        xd = 10;
        yd = 2*w;
        xd_dw = 0;
        yd_dw = 2;
        xd_ddw = 0; yd_ddw = 0;
    case 2
        %---------------正弦曲线路径--------------
        xd = 8*cos(0.2*w) + 2*w;
        yd = 2*w;
        xd_dw = -1.6*sin(0.2*w)+2;
        yd_dw = 2;
        xd_ddw = -0.32*cos(0.2*w); yd_ddw = 0;
    case 3
        %---------------圆形路径------------------
        R = 10;
        xd = R*cos(w)+11; 
        yd = R*sin(w)+11; 
        xd_dw = -10*sin(w); 
        yd_dw = 10*cos(w);
        xd_ddw = -10*cos(w); yd_ddw = -10*sin(w);
    case 4
        %---------------混合路径------------------
        xd = 10;
        yd = 2*w;
        xd_dw = 0; xd_ddw = 0;
        yd_dw = 2; yd_ddw = 0;
        if w >= 5 
            R = 5;
            xd = R*cos(w-5)+10-R;
            yd = R*sin(w-5)+10;
            xd_dw = -R*sin(w-5); xd_ddw = -R*cos(w-5);
            yd_dw = R*cos(w-5); yd_ddw = -R*sin(w-5);
        end
        if w >= 5+2*pi
            xd = 10;
            yd = 2*(w-5-2*pi)+10;
            xd_dw = 0; xd_ddw = 0;
            yd_dw = 2; yd_ddw = 0;
        end     
    case 5
        xd = 4*cos(w);
        yd = 3*w;
        xd_dw = -4*sin(w);
        yd_dw = 3;
        xd_ddw = -4*cos(w); yd_ddw = 0;
    case 6
        xd = w;
        yd = w;
        xd_dw = 1;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
        if w>=15
            xd = 30-w;
            yd = w;
            xd_dw = -1;
            yd_dw = 1;
            xd_ddw = 0; yd_ddw = 0;
        end
    case 7
        xd = w;
        yd = w;
        xd_dw = 1;
        yd_dw = 1;
        xd_ddw = 0; yd_ddw = 0;
end

%---------------------------------------------------------
kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; % 曲线路径的曲率

psif = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
  
xe = cos(psif) * (x - xd) + sin(psif) * (y - yd);  % USV与虚拟参考点的纵向误差
ye = -sin(psif) * (x - xd) + cos(psif) * (y - yd)-3;  % USV与虚拟参考点的横向误差

Xsf = psi + beta - psif; 
up = U*cos(Xsf) + k1*xe;   % 路径虚拟参考点的切向速度
wdot =  up / sqrt(xd_dw^2 + yd_dw^2);  % 路径参数的更新
w = euler2(wdot, w, ts);

Xsf_d = -atan2(ye, delta); 

end

