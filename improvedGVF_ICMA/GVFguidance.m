function [wd,e,xd,yd,y,e_alphap,Xd_dot]  = GVFguidance( pc,Ud,Vchat,Vchat_dot,p,p_dot,X,ts )
%GVFGUIDANCE wd  = GVFguidance( Ud,Vchat,Vchat_dot,p,p_dot,X,ts )
%   向量场制导，返回期望转艏率
% 海流以及侧滑速度均已知
persistent theta k1 k2 kr L D
if isempty(theta)
   theta = 0; 
   k1 = 1; k2 = 1;
   kr = 1;
   L = 0.01; D = 1;
end
% 常用变量计算
mn = [cos(X) sin(X)]'; mc = Ud*mn+Vchat;
mcn = mc/norm(mc); 
x = p(1); y = p(2); % 无人艇位置
E = [0 -1; 1 0];
K = diag([k1 k2]);
% 真实路径信息
switch pc
    case 1
       % 圆形路径
       R = 6; O = 10;
       f1 = R*cos(theta)+O;
       f2 = R*sin(theta)+O;
       xd = f1; yd = f2;
       df1_theta = -R*sin(theta);
       df2_theta =  R*cos(theta);
       ddf1_theta = -R*cos(theta);
       ddf2_theta = -R*sin(theta);
    case 2
       % 正弦路径
       f1 = theta;
       f2 = 20*sin(pi*theta/25);
       xd = f1; yd = f2;
       df1_theta = 1;
       df2_theta =  20*cos(pi*theta/25)*pi/25;
       ddf1_theta = 0;
       ddf2_theta = -20*sin(theta)*(pi/25)^2;
    case 3
       % 三叶结(trefoil knot)
       a1 = 0.02; b1 = 10; c1 = 0.03; d1 = 20; e1 = 0;
       a2 = 0.02; b2 = 10; c2 = 0.03; d2 = 20; e2 = 0;
       f1 = cos(a1*theta)*(b1*cos(c1*theta)+d1)+e1;
       f2 = sin(a2*theta)*(b2*cos(c2*theta)+d2)+e2;
       xd = f1; yd = f2;
       df1_theta = a1*sin(a1*theta)*(b1*cos(c1*theta)+d1)+...
                   cos(a1*theta)*(-b1*c1*sin(c1*theta));
       df2_theta = a2*cos(a2*theta)*(b2*cos(c2*theta)+d2)+...
                   sin(a2*theta)*(-b2*c2*sin(c2*theta));
       ddf1_theta = a1^2*cos(a1*theta)*(b1*cos(c1*theta)+d1)+...
                    a1*sin(a1*theta)*(-b1*c1*sin(c1*theta))-...
                    -a1*sin(a1*theta)*(-b1*c1*sin(c1*theta))+...
                    cos(a1*theta)*(-b1*c1^2*cos(c1*theta));
       ddf2_theta = -a2^2*sin(a2*theta)*(b2*cos(c2*theta)+d2)+...
                    a2*cos(a2*theta)*(-b2*c2*sin(c2*theta))+...
                    a2*cos(a2*theta)*(-b2*c2*sin(c2*theta))+...
                    sin(a2*theta)*(-b2*c2^2*cos(c2*theta));
    case 4
       % Lissajous曲线
       a1 = 20; b1 = 5; c1 = 0.1; d1 = 0;
       a2 = 20; b2 = 7; c2 = 0.2; d2 = 0;
       f1 = a1*cos(b1* theta+c1)+d1;
       f2 = a2*cos(b2*theta+c2)+d2;
       xd = f1; yd = f2;
       df1_theta = -a1*b1*sin(b1*theta+c1);
       df2_theta = -a2*b2*sin(b2*theta+c2);
       ddf1_theta = -a1*b1^2*cos(b1*theta+c1);
       ddf2_theta = -a2*b2^2*cos(b2*theta+c2);  
    case 5
       % 直线
       f1 = theta;
       f2 = theta;
       xd = f1; yd = f2;
       df1_theta = 1;
       df2_theta = 1;
       ddf1_theta = 0;
       ddf2_theta = 0;  
end
% 延伸至三维路径(虚拟路径)
phi1 = L*(x-f1); phi2 = L*(y-f2);

% 制导向量场
n1 = L*[1 0 -df1_theta]';
n2 = L*[0 1 -df2_theta]';
N = [n1 n2];
taun = L^2*[df1_theta df2_theta 1]';
e = [phi1 phi2]';
md = D*taun-N*K*e;
mpd = [md(1) md(2)]';
mdn = md/norm(mpd);
mpdn = mpd/norm(mpd);
% 制导率
theta_dot = norm(mcn)*mdn(3);
theta = euler2(theta_dot,theta,ts);
P_dot = [p_dot;theta_dot];
mpd_dot = L^2*[-k1 0 D*ddf1_theta+k1*df1_theta;
                0 -k2 D*ddf2_theta+k2*df1_theta]*P_dot;
mpdn_dot = -E*(mpdn*mpdn')*E*mpd_dot/norm(mpd);
Xd_dot = -mpdn'*E*mpdn_dot;

wd1 = Xd_dot-kr*mcn'*E*mpdn;
wd = (norm(mc)*wd1-mcn'*E*Vchat_dot)/(Ud*mcn'*mn);

e_alphap = asin(mcn'*E*mpdn);

y = theta;

end

