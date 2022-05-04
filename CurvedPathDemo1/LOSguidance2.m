function [sys,x0,str,ts,simStateCompliance] = LOSguidance2(t,x,u,flag)
switch flag
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
  case 1
    sys=mdlDerivatives(t,x,u);
  case {2,4,9}
    sys=[];
  case 3
    sys=mdlOutputs(t,x,u);
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
% end sfuntmpl
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 1;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 3;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = 0;

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
w = x; psi = u(6); uu = u(1); v = u(2); xx = u(4); y = u(5);
%---------------直线路径----------------------------------
K1 = 2;  % 路径虚拟参考点速度更新参数
xd = 30;
yd = 5*w;
xd_dw = 0;
yd_dw = 5;
% xd_ddw = 0;
% yd_ddw = 0;
% kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3;
psif = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
xxsf = psif-psi;
s = cos(psif)*(xx-xd)+sin(psif)*(y-yd);  % USV与虚拟参考点的纵向误差
Ud = sqrt(uu^2+v^2)*cos(xxsf)+K1*s;   % 路径虚拟参考点的切向速度
wdot =  Ud/sqrt(xd_dw^2+yd_dw^2);  % 路径参数的更新
xdot = wdot;
sys=xdot;
% end mdlDerivatives

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
uu = u(1); v = u(2); r = u(3); xx = u(4); y = u(5); psi = u(6);
w = x;
%---------------直线路径----------------------------------
xd = 30;
yd = 5*w;
xd_dw = 0;
yd_dw = 5;
% xd_ddw = 0;
% yd_ddw = 0;
% kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3;
psif = atan2(yd_dw,xd_dw);  % 路径虚拟参考点切线与x轴夹角
e = -sin(psif)*(xx-xd)+cos(psif)*(y-yd);  % USV与虚拟参考点的横向误差
deta = 3;
xxd = atan2(e,deta);  
psid = psif - xxd;
sys(1) = psid;
sys(2) = xd;
sys(3) = yd;
% end mdlOutputs

