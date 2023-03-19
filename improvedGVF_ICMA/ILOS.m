function [ psid,  ye ] = ILOS( delta, x, y, psi, u, ts)
%ILOS1 [ psid, afak, ye, yeint, pointer ] = ILOS1( point_database, R, x, y, yeint, pointer )
%   calculate disaired course angle according to ALOS
% reference: <<Integral LOS Control for Path Following of Underactuated Marine
% Surface Vessels in the Presence of Constant Ocean Currents>> by Even B?rhaug, A. Pavlov and Kristin Y. Pettersen
persistent w
if isempty(w)
    w = 0;
end
% curved path
switch pc
    case 1
        %---------------直线路径-----------------
        xd = 10;
        yd = 5*w;
        xd_dw = 0;
        yd_dw = 5;
        xd_ddw = 0; yd_ddw = 0;
    case 2
        %---------------正弦曲线路径--------------
        xd = 8*cos(0.2*w) + 2*w;
        yd = 5*w;
        xd_dw = -1.6*sin(0.2*w)+2;
        yd_dw = 5;
        xd_ddw = -0.32*cos(0.2*w); yd_ddw = 0;
    case 3
        %---------------圆形路径------------------
        R = 5;
        xd = R*cos(w)+10; 
        yd = R*sin(w)+7; 
        xd_dw = -R*sin(w); 
        yd_dw = R*cos(w);
        xd_ddw = -R*cos(w); yd_ddw = -R*sin(w);
    case 4
        %---------------混合路径------------------
        xd = 15;
        yd = 5*w;
        xd_dw = 0;
        yd_dw = 5;
        xd_ddw = 0; yd_ddw = 0;
        if w >= 6
            R = 5;
            xd = R*cos(w-6)+15-R;
            yd = R*sin(w-6)+30;
            xd_dw = -R*sin(w-6);
            yd_dw = R*cos(w-6);
            xd_ddw = -R*cos(w-6); yd_ddw = -R*sin(w-6);
        end
        if w >= 6+pi
            xd = 5;
            yd = 5*(-w+6+pi)+30;
            xd_dw = 0;
            yd_dw = -5;
            xd_ddw = 0; yd_ddw = 0;
        end        
        if w >= 6+6+pi
            R = 5;
            xd = R*cos(w-12)+5+R;
            yd = R*sin(w-12);
            xd_dw = -R*sin(w-12);
            yd_dw = R*cos(w-12);
            xd_ddw = -R*cos(w-12); yd_ddw = -R*sin(w-12);
        end  
        if w >= 6+6+pi+pi
            xd = 15;
            yd = 5*(w-12-2*pi);
            xd_dw = 0;
            yd_dw = 5;
            xd_ddw = 0; yd_ddw = 0;
        end  
    case 5
        %---------------正弦曲线路径--------------
        xd = 5*cos(0.8*w);
        yd = 5*w;
        xd_dw = -4*sin(0.8*w);
        yd_dw = 5;
        xd_ddw = -3.2*cos(0.8*w); yd_ddw = 0;
end

%---------------------------------------------------------
kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3; % 曲线路径的曲率
psif = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
  
xe = cos(psif) * (x - xd) + sin(psif) * (y - yd);  % USV与虚拟参考点的纵向误差
ye = -sin(psif) * (x - xd) + cos(psif) * (y - yd);  % USV与虚拟参考点的横向误差

persistent yeint sigma
if isempty(yeint)
    yeint = ye;
    sigma = 0.01;
end

% calculate adaptive integral term
yeintdot = sigma*delta*ye/((ye+sigma*yeint)^2+delta^2);
yeint = euler2(yeintdot,yeint,ts);
% calculate disaired course angle according LOS law
psid = psif-atan2((ye+sigma*yeint)/delta,1);
% update law
up = u*cos(psi-psif)+k*xe;
w_dot = up/sqrt(xd_dw^2+yd_dw^2);
w = euler2(w_dot,w,ts);
end

